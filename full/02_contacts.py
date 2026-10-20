"""Full run, step 2: ligand-water contacts + pocket data, streamed with NO structure files written to disk.

Each worker downloads <PDB>.cif.gz into memory, parses it, extracts what later steps need, and drops it.
Output, per chunk of entries (resumable: finished chunks are skipped on restart), in ./data/contacts/:
  contacts_NNNN.csv.gz  one row per (ligand instance, water within 5 A): interaction type, distances, water B, partners
  pockets_NNNN.pkl.gz   per ligand instance: ligand atoms, waters within 8 A, and the CA trace of the contacted protein
                        chain (used by the later superposition/site step, so the structures never need re-downloading)
  failures_NNNN.csv     entries that could not be processed

Usage (from this folder):
  python 02_contacts.py --limit 40                 # smoke test
  python 02_contacts.py --top-targets 20           # dry run on the 20 most populated UniProts
  python 02_contacts.py                            # everything
"""
import os, sys, gzip, time, pickle, argparse, requests, gemmi, numpy as np, pandas as pd
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data")
DIRECT, PROX, METAL_D, SHELL = 3.5, 5.0, 2.8, 8.0
BB = {"N", "CA", "C", "O", "OXT"}
METAL_RES = {"MG", "ZN", "MN", "CA", "FE", "FE2", "CO", "NI", "CU", "CD", "NA", "K", "HG", "CU1", "CU2", "MO", "W"}
MAX_INST = 6        # ligand copies kept per (entry, ligand) in the pocket file
_session = None

def fetch(pdb, tries=5):
    global _session
    if _session is None: _session = requests.Session()
    for i in range(tries):
        try:
            r = _session.get(f"https://files.rcsb.org/download/{pdb}.cif.gz", timeout=180); r.raise_for_status()
            return gzip.decompress(r.content).decode()       # in memory only
        except Exception as e:
            err = e; time.sleep(2 * (i + 1))
    raise err

def is_aa(name):
    ri = gemmi.find_tabulated_residue(name); return bool(ri and ri.is_amino_acid())

def analyse(pdb, ligs, ent2up):
    """ligs: set of ligand comp ids; ent2up: {entity_id: (uniprot, organism)}"""
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(fetch(pdb)).sole_block())
    st.setup_entities(); st.remove_hydrogens()
    model = st[0]; cell = st.cell; ns = gemmi.NeighborSearch(model, cell, SHELL + 1).populate()
    prot = [a.b_iso for ch in model for r in ch if is_aa(r.name) for a in r]
    bprot = float(np.mean(prot)) if prot else np.nan
    rows, pockets, ninst = [], [], Counter()
    for ch in model:
        for res in ch:
            if res.name not in ligs or ninst[res.name] >= MAX_INST: continue
            wat, cnt = {}, Counter()
            for a in res:
                for m in ns.find_atoms(a.pos, '\0', radius=SHELL):
                    cra = m.to_cra(model); rn = cra.residue.name
                    im = cell.find_nearest_image(a.pos, m.pos, gemmi.Asu.Any); d = im.dist()
                    if rn == "HOH":
                        if d > SHELL: continue
                        k = (cra.chain.name, cra.residue.seqid.num, cra.atom.altloc)
                        w = wat.setdefault(k, dict(cra=cra, pos=cell.find_nearest_pbc_position(a.pos, cra.atom.pos, m.image_idx), dmin=99., dpol=99., sym=im.sym_idx != 0))
                        w["dmin"] = min(w["dmin"], d)
                        if a.element.name in ("N", "O"): w["dpol"] = min(w["dpol"], d)
                    elif m.image_idx == 0 and d <= 5.0 and is_aa(rn): cnt[cra.chain.name] += 1
            if not cnt: continue                             # no protein contact in the asymmetric unit
            fit = model[cnt.most_common(1)[0][0]].get_polymer()
            ent = st.get_entity_of(fit); eid = ent.name if ent else ""
            up, org = ent2up.get(eid, ("", ""))
            ninst[res.name] += 1
            for k, w in wat.items():
                if w["dmin"] > PROX: continue
                wp = w["pos"]; pd_, pname, pbb, md, mname = 99., "", None, 99., ""
                for m in ns.find_atoms(wp, '\0', radius=4.0):
                    cra = m.to_cra(model); rn = cra.residue.name
                    if rn == "HOH" or (cra.chain.name == ch.name and cra.residue.seqid == res.seqid and rn == res.name): continue
                    d = cell.find_nearest_image(wp, m.pos, gemmi.Asu.Any).dist()
                    if is_aa(rn) and cra.atom.element.name in ("N", "O") and d < pd_:
                        pd_, pname, pbb = d, f"{rn}{cra.residue.seqid.num}:{cra.atom.name}", cra.atom.name in BB
                    elif cra.atom.element.is_metal or rn in METAL_RES:
                        if d < md: md, mname = d, f"{rn}{cra.residue.seqid.num}"
                direct = w["dpol"] <= DIRECT; c = w["cra"]
                itype = ("direct_water_metal" if direct and md <= METAL_D else "direct_water_protein" if direct and pd_ <= DIRECT
                         else "direct_water_only" if direct else "proximal_no_polar_contact")
                rows.append(dict(pdb=pdb, ligand=res.name, lig_chain=ch.name, lig_seqid=res.seqid.num, uniprot=up, organism=org, entity_id=eid, fit_chain=cnt.most_common(1)[0][0],
                    water=f"{c.chain.name}{c.residue.seqid.num}", water_occ=c.atom.occ, water_B=c.atom.b_iso, protein_meanB=bprot,
                    B_ratio=c.atom.b_iso / bprot if bprot else np.nan, d_lig_min=round(w["dmin"], 2),
                    d_lig_polar=round(w["dpol"], 2) if w["dpol"] < 99 else np.nan, interaction=itype, protein_partner=pname,
                    d_protein=round(pd_, 2) if pd_ < 99 else np.nan, partner_is_backbone=pbb, metal=mname,
                    d_metal=round(md, 2) if md < 99 else np.nan, symmetry_water=w["sym"]))
            ca = [(r.seqid.num, gemmi.find_tabulated_residue(r.name).one_letter_code.upper() if gemmi.find_tabulated_residue(r.name) else "X", r["CA"][0].pos)
                  for r in fit if r.find_atom("CA", '*')]
            pockets.append(dict(pdb=pdb, ligand=res.name, lig_chain=ch.name, lig_seqid=res.seqid.num, uniprot=up, organism=org, entity_id=eid,
                protein_meanB=bprot, lig=[(a.element.name, (a.pos.x, a.pos.y, a.pos.z)) for a in res],
                wat=[((w["pos"].x, w["pos"].y, w["pos"].z), w["cra"].atom.b_iso, w["cra"].atom.occ, w["sym"]) for w in wat.values()],
                ca_num=np.array([c[0] for c in ca], np.int32), ca_aa="".join(c[1] for c in ca),
                ca_xyz=np.array([(c[2].x, c[2].y, c[2].z) for c in ca], np.float32)))
    return rows, pockets

def job(args):
    pdb, ligs, ent2up = args
    try: return pdb, *analyse(pdb, set(ligs), ent2up), None
    except Exception as e: return pdb, [], [], f"{type(e).__name__}: {e}"

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12); ap.add_argument("--chunk", type=int, default=500)
    ap.add_argument("--top-targets", type=int, default=0, help="only entries containing one of the N most populated human UniProts")
    ap.add_argument("--limit", type=int, default=0, help="only the first N entries (smoke test)")
    ap.add_argument("--out", default=os.path.join(DATA, "contacts"))
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    L = pd.read_csv(f"{DATA}/entries_ligands.csv"); E = pd.read_csv(f"{DATA}/entity_map.csv", dtype=str).fillna("")
    if a.top_targets:
        top = set(pd.read_csv(f"{DATA}/target_counts.csv").uniprot.head(a.top_targets))
        L = L[L.pdb.isin(set(E[E.uniprot.isin(top)].pdb))]
    ids = sorted(L.pdb.unique());
    if a.limit: ids = ids[:a.limit]
    lig = L.groupby("pdb").ligand.apply(lambda s: sorted(set(s))).to_dict()
    ent = {p: {r.entity_id: (r.uniprot, r.organism) for r in g.itertuples()} for p, g in E[E.pdb.isin(ids)].groupby("pdb")}
    chunks = [ids[i:i + a.chunk] for i in range(0, len(ids), a.chunk)]
    print(f"{len(ids)} entries in {len(chunks)} chunks, {a.workers} workers; output {a.out}", flush=True)
    t0 = time.time(); done = 0
    with ProcessPoolExecutor(a.workers) as ex:
        for ci, chunk in enumerate(chunks):
            tag = f"{ci:04d}"
            if os.path.exists(f"{a.out}/contacts_{tag}.csv.gz") and os.path.exists(f"{a.out}/pockets_{tag}.pkl.gz"):
                done += len(chunk); print(f"chunk {tag} already done, skipping", flush=True); continue
            rows, pockets, fails = [], [], []
            futs = [ex.submit(job, (p, lig[p], ent.get(p, {}))) for p in chunk]
            for f in as_completed(futs):
                pdb, r, pk, err = f.result(); rows += r; pockets += pk
                if err: fails.append(dict(pdb=pdb, error=err))
            pd.DataFrame(rows).to_csv(f"{a.out}/contacts_{tag}.csv.gz.tmp", index=False, compression="gzip")
            with gzip.open(f"{a.out}/pockets_{tag}.pkl.gz.tmp", "wb") as fh: pickle.dump(pockets, fh, protocol=4)
            if fails: pd.DataFrame(fails).to_csv(f"{a.out}/failures_{tag}.csv", index=False)
            os.replace(f"{a.out}/pockets_{tag}.pkl.gz.tmp", f"{a.out}/pockets_{tag}.pkl.gz")
            os.replace(f"{a.out}/contacts_{tag}.csv.gz.tmp", f"{a.out}/contacts_{tag}.csv.gz")   # written last = chunk complete
            done += len(chunk); el = time.time() - t0
            print(f"chunk {tag}: {len(rows)} pairs, {len(pockets)} pockets, {len(fails)} failures | {done}/{len(ids)} entries, {el/60:.1f} min elapsed, ETA {el/done*(len(ids)-done)/60:.0f} min" if done else "", flush=True)
    print("finished", flush=True)
