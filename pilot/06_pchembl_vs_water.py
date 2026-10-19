"""Pilot step 10: pChEMBL vs ligand-water interaction metrics (per target, Spearman; size-corrected via residuals on heavy atoms)."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
w = pd.read_csv("water_contacts.csv"); SS = pd.read_csv("site_structure_table.csv"); S = pd.read_csv("water_sites.csv")
P = pd.read_csv("ligand_pchembl.csv")
cons = set(S[S.frac_present >= 0.5].site)
g = w[w.symmetry_water == False].groupby(["target", "pdb", "ligand"])
M = g.agg(n_direct_water=("interaction", lambda s: s.str.startswith("direct").sum()),
          n_bridging=("interaction", lambda s: (s == "direct_water_protein").sum()),
          mean_Bratio_direct=("B_ratio", "mean"), mw=("mw", "first"), heavy_atoms=("heavy_atoms", "first"), resolution=("resolution", "first")).reset_index()
disp = SS[(SS.status == "displaced_by_ligand") & SS.site.isin(cons)].groupby(["pdb", "ligand"]).size().rename("n_displaced_conserved")
M = M.merge(disp, on=["pdb", "ligand"], how="left").fillna({"n_displaced_conserved": 0})
M = M.merge(P[["target", "ligand", "pchembl", "n_act"]], on=["target", "ligand"])
# one row per ligand: median across its structures
L = M.groupby(["target", "ligand"]).median(numeric_only=True).reset_index()
M.to_csv("pdb_ligand_pchembl_water.csv", index=False); L.to_csv("ligand_level_pchembl_water.csv", index=False)
print("ligands with pChEMBL + water data:", L.groupby("target").size().to_dict(), "| structures:", M.groupby("target").size().to_dict())
def resid(y, x):
    ok = ~(np.isnan(y) | np.isnan(x)); c = np.polyfit(x[ok], y[ok], 1); return y - np.polyval(c, x)
print(f"\n{'target':10s} {'metric':24s} {'rho':>6s} {'p':>8s} {'rho(size-corr)':>14s} {'p':>8s}   n")
for t, d in list(L.groupby("target")) + [("ALL", L)]:
    for m in ["n_direct_water", "n_bridging", "n_displaced_conserved", "mean_Bratio_direct"]:
        x = d[m]; ok = x.notna() & d.pchembl.notna()
        if ok.sum() < 8 or x[ok].nunique() < 2: continue
        r, p = spearmanr(x[ok], d.pchembl[ok])
        pc = resid(d.pchembl.values.astype(float), d.heavy_atoms.values.astype(float))
        r2, p2 = spearmanr(x[ok], pc[ok.values])
        print(f"{t:10s} {m:24s} {r:6.2f} {p:8.1e} {r2:14.2f} {p2:8.1e}  {ok.sum()}")
r, p = spearmanr(L.heavy_atoms, L.pchembl); print("\nreference: pChEMBL vs heavy atoms rho=%.2f p=%.1e" % (r, p))
print("\nmedian pChEMBL by displaced conserved water (0 vs >=1):"); print(L.assign(disp=L.n_displaced_conserved >= 1).groupby(["target", "disp"]).pchembl.agg(["median", "size"]).round(2))
