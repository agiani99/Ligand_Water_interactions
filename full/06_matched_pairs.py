"""Full run, step 6: matched-pair test of water displacement vs potency.

Unit = a pair of DIFFERENT ligands of the same UniProt at the same conserved water site, where
  ligand A displaces the water (displaced_by_ligand) and ligand B keeps it (present),
  A and B are structurally similar (Morgan-2 Tanimoto >= --tc), crystallised at similar resolution (|dres| <= --dres),
  and both have a pChEMBL value for that target.
Effect = pChEMBL(A) - pChEMBL(B)  ( >0 means the water-displacing ligand is more potent ).
Because displacers are usually bigger, the effect is also reported after adjusting for the heavy-atom difference.
Inference is by cluster bootstrap over UniProts (pairs within a target are not independent).
"""
import os, argparse, warnings, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import wilcoxon
RDLogger.DisableLog("rdApp.*"); warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data"); R = os.path.join(DATA, "results"); os.makedirs(R, exist_ok=True)
ap = argparse.ArgumentParser(); ap.add_argument("--tc", type=float, default=0.6); ap.add_argument("--dres", type=float, default=0.5)
ap.add_argument("--min-occ", type=float, default=0.5); ap.add_argument("--min-obs", type=int, default=10); a = ap.parse_args()

S = pd.read_csv(f"{DATA}/sites/sites.csv"); O = pd.read_csv(f"{DATA}/sites/site_structure.csv.gz")
P = pd.read_csv(f"{DATA}/chembl/ligand_pchembl.csv"); L = pd.read_csv(f"{DATA}/entries_ligands.csv").drop_duplicates("ligand")[["ligand", "smiles"]]
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048); FP = {}
for lig, smi in L.values:
    m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if m is None and isinstance(smi, str): m = Chem.MolFromSmiles(smi, sanitize=False);
    try:
        if m is not None:
            if not m.GetNumAtoms(): continue
            m.UpdatePropertyCache(strict=False); Chem.GetSymmSSSR(m); FP[lig] = gen.GetFingerprint(m)
    except Exception: pass

cons = S[(S.frac_present >= a.min_occ) & (S.n_obs >= a.min_obs)][["site", "frac_present"]]
O = O[O.status.isin(["present", "displaced_by_ligand"]) & O.site.isin(cons.site)].merge(cons, on="site").merge(P[["uniprot", "ligand", "pchembl"]], on=["uniprot", "ligand"])
O = O[O.ligand.isin(FP)]
print(f"conserved-site observations with pChEMBL: {len(O):,} ({O.uniprot.nunique()} UniProts, {O.site.nunique()} sites)")

rows = []
for site, g in O.groupby("site"):
    A = g[g.status == "displaced_by_ligand"]; B = g[g.status == "present"]
    both = set(A.ligand) & set(B.ligand)                      # ligand seen on both sides in this site: ambiguous, drop
    A, B = A[~A.ligand.isin(both)], B[~B.ligand.isin(both)]
    if A.empty or B.empty: continue
    Aa = A.groupby("ligand").agg(res=("resolution", "median"), hac=("heavy_atoms", "first"), p=("pchembl", "first")); Bb = B.groupby("ligand").agg(res=("resolution", "median"), hac=("heavy_atoms", "first"), p=("pchembl", "first"))
    for la, ra in Aa.iterrows():
        sims = DataStructs.BulkTanimotoSimilarity(FP[la], [FP[lb] for lb in Bb.index])
        for lb, tc in zip(Bb.index, sims):
            rb = Bb.loc[lb]
            if tc >= min(a.tc, 0.5) and abs(ra.res - rb.res) <= max(a.dres, 1.0):
                rows.append(dict(uniprot=g.uniprot.iloc[0], site=site, frac_present=g.frac_present.iloc[0], lig_displace=la, lig_keep=lb, tanimoto=tc, res_displace=ra.res, res_keep=rb.res,
                                 dres=ra.res - rb.res, hac_displace=ra.hac, hac_keep=rb.hac, dhac=ra.hac - rb.hac, p_displace=ra.p, p_keep=rb.p, dpchembl=ra.p - rb.p))
Pr = pd.DataFrame(rows)
if Pr.empty: raise SystemExit("no matched pairs at these thresholds")
Pr = Pr.sort_values(["frac_present", "tanimoto"], ascending=False).drop_duplicates(["lig_displace", "lig_keep", "uniprot"])      # a ligand pair counts once per target
Pr.to_csv(f"{R}/matched_pairs_all_candidates.csv", index=False)

def estimate(d, boot=2000, seed=0):
    rng = np.random.default_rng(seed); groups = [g for _, g in d.groupby("uniprot")]; out = {}
    def stats(parts):
        x = pd.concat(parts); m = x.dpchembl.mean()
        # intercept (uncentred covariates) = effect at equal ligand size and equal resolution
        X = np.column_stack([np.ones(len(x)), x.dhac, x.dres]); b = np.linalg.lstsq(X, x.dpchembl.values, rcond=None)[0]
        return m, b[0], b[1]
    m, adj, slope = stats(groups); bs = np.array([stats([groups[i] for i in rng.integers(0, len(groups), len(groups))]) for _ in range(boot)])
    ci = np.percentile(bs, [2.5, 97.5], axis=0); up = d.groupby("uniprot").dpchembl.mean()
    return dict(n_pairs=len(d), n_ligand_pairs=len(d.drop_duplicates(["lig_displace", "lig_keep"])), n_uniprots=len(groups), mean_dpchembl=m, ci_lo=ci[0, 0], ci_hi=ci[1, 0],
                adj_dpchembl=adj, adj_ci_lo=ci[0, 1], adj_ci_hi=ci[1, 1], per_hac_slope=slope, frac_positive=(d.dpchembl > 0).mean(), mean_dhac=d.dhac.mean(),
                uniprot_wilcoxon_p=wilcoxon(up).pvalue if len(up) > 5 and (up != 0).any() else np.nan)

# main estimate + sensitivity grid
def build(tc, dres):
    d = Pr[(Pr.tanimoto >= tc) & (Pr.dres.abs() <= dres)]; return d
res = []
for tc in [a.tc, 0.5, 0.6, 0.7, 0.8]:
    for dres in [0.3, 0.5, 1.0]:
        d = build(tc, dres)
        if d.uniprot.nunique() >= 5: res.append(dict(tc=tc, dres=dres, **estimate(d, boot=500)))
for dh in [None, 3, 1]:
    d = build(a.tc, a.dres); d = d if dh is None else d[d.dhac.abs() <= dh]
    if d.uniprot.nunique() >= 5: res.append(dict(tc=a.tc, dres=a.dres, max_abs_dhac=dh, **estimate(d, boot=1000)))
Sg = pd.DataFrame(res).drop_duplicates().round(3); Sg.to_csv(f"{R}/matched_pairs_sensitivity.csv", index=False)
main = estimate(build(a.tc, a.dres), boot=5000)
print(f"\nMAIN (Tanimoto >= {a.tc}, |dres| <= {a.dres}, sites occupied >= {a.min_occ:.0%}):")
for k, v in main.items(): print(f"  {k}: {round(v, 3) if isinstance(v, float) else v}")
print("\nSensitivity:\n" + Sg[["tc", "dres", "max_abs_dhac", "n_pairs", "n_uniprots", "mean_dpchembl", "ci_lo", "ci_hi", "adj_dpchembl", "adj_ci_lo", "adj_ci_hi", "mean_dhac"]].to_string(index=False)) if "max_abs_dhac" in Sg else print(Sg.to_string(index=False))

d = build(a.tc, a.dres); fig, ax = plt.subplots(1, 3, figsize=(15, 4))
ax[0].hist(d.dpchembl, bins=40); ax[0].axvline(0, c="k", lw=.8); ax[0].axvline(d.dpchembl.mean(), c="r"); ax[0].set(xlabel="pChEMBL(displacing) - pChEMBL(retaining)", ylabel="ligand pairs", title=f"{len(d)} matched pairs, {d.uniprot.nunique()} targets")
ax[1].scatter(d.dhac, d.dpchembl, s=8, alpha=.4); ax[1].axhline(0, c="k", lw=.5); ax[1].axvline(0, c="k", lw=.5); ax[1].set(xlabel="heavy atoms (displacing - retaining)", ylabel="delta pChEMBL", title="Size confounding")
up = d.groupby("uniprot").agg(m=("dpchembl", "mean"), n=("dpchembl", "size")).sort_values("n", ascending=False).head(25)
ax[2].barh(range(len(up)), up.m); ax[2].set(yticks=range(len(up)), yticklabels=[f"{u} (n={n})" for u, n in zip(up.index, up.n)], xlabel="mean delta pChEMBL", title="Top targets by number of pairs"); ax[2].invert_yaxis(); ax[2].axvline(0, c="k", lw=.5)
plt.tight_layout(); plt.savefig(f"{R}/matched_pairs.png", dpi=150)
