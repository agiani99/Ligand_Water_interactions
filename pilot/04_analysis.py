"""Pilot step 4: water-site clustering, displacement calls, stats and plots."""
import pickle, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, mannwhitneyu

SITE_R, DISP_R, DIRECT = 1.2, 2.5, 3.5
obs = pickle.load(open("obs.pkl", "rb"))
ent = pd.read_csv("entries_ligands.csv").drop_duplicates(["pdb", "ligand"])
cont = pd.read_csv("water_contacts.csv")
bprot = cont.drop_duplicates("pdb").set_index("pdb").protein_meanB.to_dict()
meta = ent.set_index(["pdb", "ligand"])

site_rows, ss_rows = [], []
for target, lst in obs.items():
    pts, owner, bs, oc = [], [], [], []
    for i, o in enumerate(lst):
        for q, b, o_ in o["wat"]: pts.append(q); owner.append(i); bs.append(b); oc.append(o_)
    pts, owner, bs = np.array(pts), np.array(owner), np.array(bs)
    tree = cKDTree(pts); free = np.ones(len(pts), bool); sid = 0
    ligxyz = [np.array([q for e, q in o["lig"]]) for o in lst]
    ligpol = [np.array([q for e, q in o["lig"] if e in ("N", "O")]).reshape(-1, 3) for o in lst]
    ligtrees = [cKDTree(x) for x in ligxyz]
    sites = []
    while True:
        nb = [len([j for j in n if free[j]]) for n in tree.query_ball_point(pts, SITE_R)] if sid == 0 else None
        break
    # greedy densest-first clustering
    nbrs = tree.query_ball_point(pts, SITE_R)
    while True:
        cnt = np.array([sum(free[j] for j in nbrs[k]) if free[k] else 0 for k in range(len(pts))])
        k = cnt.argmax()
        if cnt[k] < 3: break
        mem = [j for j in nbrs[k] if free[j]]
        c = pts[mem].mean(0); best = {}
        for j in tree.query_ball_point(c, SITE_R):
            if free[j] and (owner[j] not in best or np.linalg.norm(pts[j]-c) < np.linalg.norm(pts[best[owner[j]]]-c)): best[owner[j]] = j
        for j in best.values(): free[j] = False
        c = pts[list(best.values())].mean(0); sites.append((c, best))
    n = len(lst)
    for s, (c, best) in enumerate(sites):
        near = [i for i in range(n) if ligtrees[i].query(c)[0] <= 5.0]
        if len(near) < 5: continue          # site not in ligand pocket for enough structures
        sname = f"{target}_W{s+1:03d}"
        npres = ndisp = 0
        for i in near:
            o = lst[i]; dl = ligtrees[i].query(c)[0]
            dpol = np.linalg.norm(ligpol[i] - c, axis=1).min() if len(ligpol[i]) else np.nan
            m = meta.loc[(o["pdb"], o["ligand"])]
            if i in best:
                j = best[i]; status = "present"; B = bs[j]; npres += 1
                dl_w = ligtrees[i].query(pts[j])[0]
                dpol_w = np.linalg.norm(ligpol[i] - pts[j], axis=1).min() if len(ligpol[i]) else np.nan
            else:
                status = "displaced_by_ligand" if dl < DISP_R else "absent"; B = np.nan; dl_w = dpol_w = np.nan
                ndisp += status == "displaced_by_ligand"
            ss_rows.append(dict(target=target, site=sname, pdb=o["pdb"], ligand=o["ligand"], resolution=m.resolution, mw=m.mw, heavy_atoms=m.heavy_atoms,
                status=status, water_B=B, B_ratio=B/bprot[o["pdb"]] if status == "present" else np.nan,
                d_ligand_site=round(dl, 2), d_ligand_water=dl_w, d_ligand_polar_water=dpol_w,
                direct_contact=bool(dpol_w <= DIRECT) if status == "present" else np.nan))
        site_rows.append(dict(target=target, site=sname, x=c[0], y=c[1], z=c[2], n_structures_pocket=len(near), n_present=npres, n_displaced=ndisp))

S = pd.DataFrame(site_rows); SS = pd.DataFrame(ss_rows)
S["frac_present"] = S.n_present / S.n_structures_pocket
S.to_csv("water_sites.csv", index=False); SS.to_csv("site_structure_table.csv", index=False)

# ---- answers
print("sites:", len(S), "| site-structure observations:", len(SS))
print(SS.groupby(["target", "status"]).size().unstack(fill_value=0))
disp = SS[SS.status == "displaced_by_ligand"]
# only count displacement for sites that are conserved (water present in >=50% of pocket structures)
cons = S[S.frac_present >= 0.5].site
dc = disp[disp.site.isin(cons)]
print("\nconserved sites (>=50% occupancy):", len(cons))
print("ligand-water displacement events at conserved sites:", len(dc), "in", dc.pdb.nunique(), "structures")
print(dc.groupby("target").agg(events=("pdb", "size"), structures=("pdb", "nunique"), sites=("site", "nunique")))
print(dc.groupby("site").size().sort_values(ascending=False).head(8))

pres = SS[SS.status == "present"].dropna(subset=["resolution"])
for col in ["water_B", "B_ratio"]:
    r, p = spearmanr(pres.resolution, pres[col]); print(f"\nSpearman resolution vs {col}: rho={r:.2f} p={p:.1e} n={len(pres)}")
frac = SS.groupby("pdb").agg(res=("resolution", "first"), pres=("status", lambda s: (s == "present").mean()), disp=("status", lambda s: (s == "displaced_by_ligand").sum()))
print("Spearman resolution vs fraction of sites with water:", spearmanr(frac.res, frac.pres))
print("Spearman resolution vs #displaced:", spearmanr(frac.res, frac.disp))

# B-factor: direct ligand contact vs no contact, resolution-corrected within site
rows = []
for s, g in pres.groupby("site"):
    a, b = g[g.direct_contact == True], g[g.direct_contact == False]
    if len(a) >= 5 and len(b) >= 5:
        coef = np.polyfit(g.resolution, g.B_ratio, 1)
        resid = g.B_ratio - np.polyval(coef, g.resolution)
        ra, rb = resid[a.index], resid[b.index]
        rows.append(dict(site=s, n_direct=len(a), n_nocontact=len(b), delta_Bratio_res_corrected=ra.median() - rb.median(), p=mannwhitneyu(ra, rb).pvalue))
D = pd.DataFrame(rows); D.to_csv("water_B_direct_vs_nocontact.csv", index=False)
print("\nB_ratio, direct ligand contact vs none (resolution-corrected), per site:"); print(D.round(3).to_string(index=False))

# ---- plots
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for t, g in pres.groupby("target"): ax[0].hist(g.water_B, bins=40, alpha=.5, label=t, density=True)
ax[0].set(xlabel="water B (Å²)", ylabel="density", title="Pocket water B-factors"); ax[0].legend()
ax[1].hist([pres[pres.direct_contact == True].B_ratio, pres[pres.direct_contact == False].B_ratio], bins=40, label=["ligand H-bond (≤3.5 Å)", "no polar contact"], density=True)
ax[1].set(xlabel="water B / mean protein B", title="Normalized B by ligand contact"); ax[1].legend()
ax[2].scatter(pres.resolution, pres.B_ratio, s=3, alpha=.25); ax[2].set(xlabel="resolution (Å)", ylabel="B ratio", title="Water B vs resolution")
plt.tight_layout(); plt.savefig("water_B_plots.png", dpi=150)
