"""Full run, step 3: per-UniProt pocket superposition, water-site clustering and presence/displacement calls.

Input : data/contacts/pockets_*.pkl.gz (from step 2), data/entries_ligands.csv
Output: data/sites/shards/<UNIPROT>_sites.csv, _obs.csv.gz, _fits.csv   (one set per UniProt, resumable)
        data/sites/sites.csv, site_structure.csv.gz, fits.csv           (combined, written when all shards exist)

Per UniProt group (>= --min-structures distinct entries):
  1. reference = best-resolution pocket among the longest chains
  2. CA atoms are matched (same residue number+aa, else Biopython global alignment), global Kabsch fit with outlier
     rejection, then a pocket-focused refit on CAs within 12 A of the ligand; pockets with pocket RMSD > 2.0 A are skipped
  3. waters within 8 A of any ligand are pooled in the reference frame and clustered greedily (1.2 A) into sites
  4. per (site, structure) with the site <= 5 A from the ligand: present | displaced_by_ligand (ligand atom < 2.5 A) | absent
Only the first ligand copy per (entry, ligand) is used.
"""
import os, sys, gzip, glob, pickle, argparse, time, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio.Align import PairwiseAligner

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data"); SH = os.path.join(DATA, "sites", "shards")
SITE_R, DISP_R, NEAR_R, DIRECT = 1.2, 2.5, 5.0, 3.5
POCKET_R, POCKET_RMSD_MAX, MIN_SITE_MEMBERS, MIN_NEAR = 12.0, 2.0, 3, 5

def kabsch(P, Q):
    """rotation/translation mapping Q onto P (both n x 3)."""
    pc, qc = P.mean(0), Q.mean(0); H = (Q - qc).T @ (P - pc)
    U, S, Vt = np.linalg.svd(H); d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, pc - R @ qc                                   # x' = R @ x + t  (apply as x @ R.T + t)

def apply(R, t, X): return np.asarray(X, float) @ R.T + t

_al = None
def match(ref, mob):
    """indices (ref_idx, mob_idx) of equivalent CA atoms."""
    rd = {n: i for i, n in enumerate(ref["ca_num"])}
    pairs = [(rd[n], j) for j, n in enumerate(mob["ca_num"]) if n in rd and ref["ca_aa"][rd[n]] == mob["ca_aa"][j]]
    if len(pairs) >= max(30, 0.6 * min(len(ref["ca_num"]), len(mob["ca_num"]))): return np.array(pairs)
    global _al
    if _al is None:
        _al = PairwiseAligner(); _al.mode = "global"; _al.match_score, _al.mismatch_score = 2, -1; _al.open_gap_score, _al.extend_gap_score = -2, -0.5
    aln = _al.align(ref["ca_aa"], mob["ca_aa"])[0]; out = []
    for (r0, r1), (m0, m1) in zip(*aln.aligned):
        out += [(r0 + k, m0 + k) for k in range(r1 - r0) if ref["ca_aa"][r0 + k] == mob["ca_aa"][m0 + k]]
    return np.array(out)

def fit(ref, mob):
    """returns (R, t, pocket_rmsd, n_pocket, reason)"""
    pr = match(ref, mob)
    if len(pr) < 20: return None, None, np.nan, 0, "too_few_matched_CA"
    P, Q = ref["ca_xyz"][pr[:, 0]].astype(float), mob["ca_xyz"][pr[:, 1]].astype(float); sel = np.ones(len(P), bool)
    for cut in (4.0, 3.0, 2.5):                              # global fit with outlier rejection
        if sel.sum() < 10: return None, None, np.nan, 0, "global_fit_failed"
        R, t = kabsch(P[sel], Q[sel]); d = np.linalg.norm(apply(R, t, Q) - P, axis=1); sel = d < cut
    lc = apply(R, t, np.array([x for _, x in mob["lig"]]).mean(0, keepdims=True))[0]
    ps = np.linalg.norm(P - lc, axis=1) < POCKET_R
    if ps.sum() < 8: return None, None, np.nan, int(ps.sum()), "too_few_pocket_CA"
    R, t = kabsch(P[ps], Q[ps]); d = np.linalg.norm(apply(R, t, Q[ps]) - P[ps], axis=1)
    keep = d < 3.0                                           # one more rejection pass inside the pocket
    if keep.sum() >= 8: R, t = kabsch(P[ps][keep], Q[ps][keep]); d = np.linalg.norm(apply(R, t, Q[ps][keep]) - P[ps][keep], axis=1)
    rmsd = float(np.sqrt((d ** 2).mean()))
    return (R, t, rmsd, int(len(d)), "ok") if rmsd <= POCKET_RMSD_MAX else (None, None, rmsd, int(len(d)), "pocket_rmsd_too_high")

def process_group(up, pockets, meta):
    ref = max([p for p in pockets if len(p["ca_num"]) >= 0.8 * max(len(q["ca_num"]) for q in pockets)], key=lambda p: -meta[(p["pdb"], p["ligand"])]["resolution"])
    fits, obs = [], []
    for p in pockets:
        if p is ref: R, t, rmsd, n, why = np.eye(3), np.zeros(3), 0.0, len(p["ca_num"]), "ok"
        else: R, t, rmsd, n, why = fit(ref, p)
        fits.append(dict(uniprot=up, pdb=p["pdb"], ligand=p["ligand"], pocket_rmsd=round(rmsd, 2) if rmsd == rmsd else np.nan, n_pocket_ca=n, status=why, is_ref=p is ref))
        if why != "ok": continue
        m = meta[(p["pdb"], p["ligand"])]; lig = apply(R, t, [x for _, x in p["lig"]]); el = np.array([e for e, _ in p["lig"]])
        wat = np.array([w[0] for w in p["wat"]]).reshape(-1, 3)
        obs.append(dict(pdb=p["pdb"], ligand=p["ligand"], res=m["resolution"], mw=m["mw"], hac=m["heavy_atoms"], bprot=p["protein_meanB"],
                        lig=lig, polar=lig[np.isin(el, ["N", "O"])], wat=apply(R, t, wat) if len(wat) else wat,
                        wB=np.array([w[1] for w in p["wat"]]), ltree=cKDTree(lig)))
    F = pd.DataFrame(fits)
    if len(obs) < 5: return up, pd.DataFrame(), pd.DataFrame(), F
    for o in obs: o["ptree"] = cKDTree(o["polar"]) if len(o["polar"]) else None
    pts = np.vstack([o["wat"] for o in obs if len(o["wat"])]); owner = np.concatenate([[i] * len(o["wat"]) for i, o in enumerate(obs)]).astype(int)
    Bw = np.concatenate([o["wB"] for o in obs if len(o["wat"])])
    if len(pts) == 0: return up, pd.DataFrame(), pd.DataFrame(), F
    tree = cKDTree(pts); nbrs = tree.query_ball_point(pts, SITE_R); cnt = np.array([len(n) for n in nbrs]); free = np.ones(len(pts), bool); sites = []
    while True:
        k = int(np.argmax(np.where(free, cnt, -1)))
        if cnt[k] < MIN_SITE_MEMBERS: break
        c = pts[[j for j in nbrs[k] if free[j]]].mean(0); best = {}
        for j in tree.query_ball_point(c, SITE_R):
            if free[j] and (owner[j] not in best or np.linalg.norm(pts[j] - c) < np.linalg.norm(pts[best[owner[j]]] - c)): best[owner[j]] = j
        if not best: best = {owner[k]: k}
        for j in best.values():
            free[j] = False
            for n in nbrs[j]: cnt[n] -= 1
        sites.append((pts[list(best.values())].mean(0), best))
    if not sites: return up, pd.DataFrame(), pd.DataFrame(), F
    C = np.array([s[0] for s in sites])
    D = np.array([o["ltree"].query(C)[0] for o in obs])      # obs x sites : nearest ligand atom to site centre
    rows, srows = [], []
    for s, (c, best) in enumerate(sites):
        near = np.where(D[:, s] <= NEAR_R)[0]
        if len(near) < MIN_NEAR: continue
        name = f"{up}_W{s + 1:04d}"; npres = ndisp = 0
        for i in near:
            o = obs[i]; dl = D[i, s]
            if i in best:
                j = best[i]; st = "present"; npres += 1; B = Bw[j]; dlw = o["ltree"].query(pts[j])[0]
                dpw = o["ptree"].query(pts[j])[0] if o["ptree"] is not None else np.nan
            else:
                st = "displaced_by_ligand" if dl < DISP_R else "absent"; B = dlw = dpw = np.nan; ndisp += st == "displaced_by_ligand"
            rows.append((name, o["pdb"], o["ligand"], o["res"], o["mw"], o["hac"], st, B, B / o["bprot"] if st == "present" else np.nan,
                         round(dl, 2), dlw, dpw, (dpw <= DIRECT) if st == "present" else np.nan))
        srows.append(dict(uniprot=up, site=name, x=c[0], y=c[1], z=c[2], ref_pdb=ref["pdb"], n_obs=len(near), n_present=npres, n_displaced=ndisp))
    cols = ["site", "pdb", "ligand", "resolution", "mw", "heavy_atoms", "status", "water_B", "B_ratio", "d_ligand_site", "d_ligand_water", "d_ligand_polar_water", "direct_contact"]
    O = pd.DataFrame(rows, columns=cols); O.insert(0, "uniprot", up)
    S = pd.DataFrame(srows); S["frac_present"] = S.n_present / S.n_obs if len(S) else []
    return up, S, O, F

def run_group(args):
    up, pockets, meta = args
    try:
        up, S, O, F = process_group(up, pockets, meta)
        S.to_csv(f"{SH}/{up}_sites.csv", index=False); O.to_csv(f"{SH}/{up}_obs.csv.gz", index=False); F.to_csv(f"{SH}/{up}_fits.csv", index=False)
        return up, len(F), int((F.status == "ok").sum()), len(S), None
    except Exception as e:
        import traceback; return up, 0, 0, 0, traceback.format_exc()[-400:]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12); ap.add_argument("--min-structures", type=int, default=5)
    ap.add_argument("--uniprots", default="", help="comma list, for validation runs"); ap.add_argument("--redo", action="store_true")
    a = ap.parse_args(); os.makedirs(SH, exist_ok=True)
    meta = pd.read_csv(f"{DATA}/entries_ligands.csv").drop_duplicates(["pdb", "ligand"]).set_index(["pdb", "ligand"]).to_dict("index")
    groups, seen = {}, set()
    for f in sorted(glob.glob(f"{DATA}/contacts/pockets_*.pkl.gz")):
        for p in pickle.load(gzip.open(f, "rb")):
            k = (p["pdb"], p["ligand"])
            if k in seen or not p["uniprot"] or "Homo sapiens" not in p["organism"] or k not in meta or len(p["ca_num"]) < 30: continue
            seen.add(k); groups.setdefault(p["uniprot"], []).append(p)
    want = set(a.uniprots.split(",")) if a.uniprots else None
    todo = {u: g for u, g in groups.items() if len({p["pdb"] for p in g}) >= a.min_structures and (want is None or u in want)
            and (a.redo or not os.path.exists(f"{SH}/{u}_fits.csv"))}
    print(f"{len(groups)} UniProt groups; {len(todo)} to process (>= {a.min_structures} structures)", flush=True)
    t0 = time.time(); n = 0
    with ProcessPoolExecutor(a.workers) as ex:
        futs = [ex.submit(run_group, (u, g, {(p["pdb"], p["ligand"]): meta[(p["pdb"], p["ligand"])] for p in g})) for u, g in sorted(todo.items(), key=lambda x: -len(x[1]))]
        for f in as_completed(futs):
            up, nf, nok, ns, err = f.result(); n += 1
            print(f"[{n}/{len(todo)}] {up}: {nok}/{nf} pockets fitted, {ns} sites" + (f"  ERROR {err}" if err else "") + f" | {(time.time() - t0) / 60:.1f} min", flush=True)
    if not want:
        out = os.path.join(DATA, "sites")
        pd.concat(pd.read_csv(f) for f in glob.glob(f"{SH}/*_sites.csv") if os.path.getsize(f) > 5).to_csv(f"{out}/sites.csv", index=False)
        pd.concat(pd.read_csv(f) for f in glob.glob(f"{SH}/*_fits.csv")).to_csv(f"{out}/fits.csv", index=False)
        first = True
        for f in glob.glob(f"{SH}/*_obs.csv.gz"):
            try: d = pd.read_csv(f)
            except pd.errors.EmptyDataError: continue
            if len(d): d.to_csv(f"{out}/site_structure.csv.gz", mode="w" if first else "a", header=first, index=False, compression="gzip"); first = False
        print("combined tables written to", out, flush=True)
