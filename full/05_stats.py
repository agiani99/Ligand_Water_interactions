"""Full run, step 5: water B-factor distributions, resolution effects, displacement counts, pChEMBL vs water interaction."""
import os, glob, warnings, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, wilcoxon, norm
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data"); R = os.path.join(DATA, "results"); os.makedirs(R, exist_ok=True)

C = pd.concat(pd.read_csv(f) for f in sorted(glob.glob(f"{DATA}/contacts/contacts_*.csv.gz")))
C = C[(C.symmetry_water == False) & C.organism.str.contains("Homo sapiens", na=False) & (C.uniprot.notna())]
L = pd.read_csv(f"{DATA}/entries_ligands.csv").drop_duplicates(["pdb", "ligand"])[["pdb", "ligand", "resolution", "mw", "heavy_atoms"]]
C = C.merge(L, on=["pdb", "ligand"], how="left")
S = pd.read_csv(f"{DATA}/sites/sites.csv"); O = pd.read_csv(f"{DATA}/sites/site_structure.csv.gz")
P = pd.read_csv(f"{DATA}/chembl/ligand_pchembl.csv")
out = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

# ---------------- 1. water B-factors
say(f"## Ligand-site waters: {len(C):,} ligand-water pairs, {C.pdb.nunique():,} entries, {C.uniprot.nunique():,} UniProts")
say("water_B percentiles (A^2):", C.water_B.quantile([.05, .25, .5, .75, .95]).round(1).to_dict())
b = C.water_B
say("Gemini bins (fraction of waters): <20: %.2f | 20-35: %.2f | 35-60: %.2f | 60-80: %.2f | >80: %.3f" %
    ((b < 20).mean(), ((b >= 20) & (b <= 35)).mean(), ((b > 35) & (b <= 60)).mean(), ((b > 60) & (b <= 80)).mean(), (b > 80).mean()))
T = C.groupby("interaction").agg(n=("water_B", "size"), median_B=("water_B", "median"), median_Bratio=("B_ratio", "median")).round(2); say(T.to_string())
T.to_csv(f"{R}/water_B_by_interaction.csv")
rb = C.groupby(pd.cut(C.resolution, [0, 1.0, 1.5, 2.0, 2.5])).agg(n=("water_B", "size"), median_B=("water_B", "median"), median_Bratio=("B_ratio", "median")).round(2); say(rb.to_string())
for col in ["water_B", "B_ratio"]:
    r, p = spearmanr(C.resolution, C[col], nan_policy="omit"); say(f"Spearman resolution vs {col}: rho={r:.2f} p={p:.1e} n={len(C):,}")
tu = C.groupby("uniprot").water_B.median(); say("median water B across UniProts: IQR %.1f-%.1f (n=%d)" % (tu.quantile(.25), tu.quantile(.75), len(tu)))

# ---------------- 2. sites, presence, displacement
cons = S[(S.frac_present >= 0.5) & (S.n_obs >= 10)]
say(f"\n## Sites: {len(S):,} sites in {S.uniprot.nunique()} UniProts; conserved (>=50% occupancy, >=10 structures): {len(cons):,} in {cons.uniprot.nunique()} UniProts")
D = O[(O.status == "displaced_by_ligand") & O.site.isin(cons.site)]
say(f"ligand-water displacement events at conserved sites: {len(D):,} in {D.pdb.nunique():,} structures, {D.uniprot.nunique()} UniProts, {D.ligand.nunique():,} distinct ligands")
say("fraction of pocket observations at conserved sites that are 'displaced':", round(float((O[O.site.isin(cons.site)].status == "displaced_by_ligand").mean()), 3))
D.groupby("uniprot").size().sort_values(ascending=False).head(15).rename("events").to_csv(f"{R}/displacement_top_uniprots.csv"); say("top UniProts by displacement events:", D.groupby("uniprot").size().sort_values(ascending=False).head(8).to_dict())
pe = O.groupby(["pdb", "ligand"]).agg(res=("resolution", "first"), pres=("status", lambda s: (s == "present").mean()), disp=("status", lambda s: (s == "displaced_by_ligand").sum())).reset_index()
r, p = spearmanr(pe.res, pe.pres); say(f"Spearman resolution vs fraction of pocket sites with a water: rho={r:.2f} p={p:.1e} (n={len(pe):,} structures)")
r, p = spearmanr(pe.res, pe.disp); say(f"Spearman resolution vs number of displaced-water calls: rho={r:.2f} p={p:.1e}")

# ---------------- 3. B ratio: direct ligand contact vs none, per site, resolution-corrected
pres = O[O.status == "present"].dropna(subset=["resolution", "B_ratio"]).copy(); rows = []
for s, g in pres.groupby("site"):
    a, n = g[g.direct_contact == True], g[g.direct_contact == False]
    if len(a) >= 5 and len(n) >= 5 and g.resolution.nunique() > 3:
        resid = g.B_ratio - np.polyval(np.polyfit(g.resolution, g.B_ratio, 1), g.resolution)
        rows.append(dict(site=s, uniprot=g.uniprot.iloc[0], n_direct=len(a), n_nocontact=len(n), delta=resid[a.index].median() - resid[n.index].median()))
Dl = pd.DataFrame(rows); Dl.to_csv(f"{R}/site_B_direct_vs_nocontact.csv", index=False)
if len(Dl):
    say(f"\n## B ratio, ligand H-bonded water vs no contact (resolution-corrected), {len(Dl)} sites: median delta {Dl.delta.median():.3f}; "
        f"{(Dl.delta < 0).mean():.0%} of sites lower for H-bonded waters; Wilcoxon p={wilcoxon(Dl.delta).pvalue:.1e}")

# ---------------- 4. pChEMBL vs water interaction
C["inst"] = C.pdb + C.ligand + C.lig_chain.astype(str) + C.lig_seqid.astype(str)
inst = C.groupby(["uniprot", "pdb", "ligand", "inst"]).agg(n_direct=("interaction", lambda s: s.str.startswith("direct").sum()),
        n_bridging=("interaction", lambda s: (s == "direct_water_protein").sum()), Bratio_direct=("B_ratio", "mean")).reset_index()
M = inst.groupby(["uniprot", "pdb", "ligand"]).mean(numeric_only=True).reset_index()
dd = D.groupby(["pdb", "ligand"]).size().rename("n_displaced_conserved")
M = M.merge(dd, on=["pdb", "ligand"], how="left").fillna({"n_displaced_conserved": 0}).merge(L[["pdb", "ligand", "heavy_atoms", "resolution"]], on=["pdb", "ligand"])
Lg = M.groupby(["uniprot", "ligand"]).median(numeric_only=True).reset_index().merge(P, on=["uniprot", "ligand"])
Lg.to_csv(f"{R}/ligand_level_pchembl_water.csv", index=False)
say(f"\n## pChEMBL vs water: {len(Lg):,} (UniProt, ligand) pairs, {Lg.uniprot.nunique()} UniProts; {int((Lg.groupby('uniprot').size() >= 15).sum())} UniProts with >=15 ligands")
res = []
for u, g in Lg.groupby("uniprot"):
    if len(g) < 15: continue
    for m in ["n_direct", "n_bridging", "n_displaced_conserved", "Bratio_direct"]:
        ok = g[m].notna() & g.pchembl.notna()
        if ok.sum() >= 15 and g[m][ok].nunique() > 1 and g.pchembl[ok].nunique() > 1:
            r, p = spearmanr(g[m][ok], g.pchembl[ok]); res.append(dict(uniprot=u, metric=m, n=int(ok.sum()), rho=r, p=p))
            h = g.heavy_atoms[ok]
            if h.notna().all() and h.nunique() > 1:                # partial Spearman controlling for ligand heavy-atom count
                rk = lambda v: v.rank().values.astype(float); x, y, z = rk(g[m][ok]), rk(g.pchembl[ok]), rk(h)
                rx, ry = x - np.polyval(np.polyfit(z, x, 1), z), y - np.polyval(np.polyfit(z, y, 1), z)
                if rx.std() > 0 and ry.std() > 0: res.append(dict(uniprot=u, metric=m + "_sizecorr", n=int(ok.sum()), rho=float(np.corrcoef(rx, ry)[0, 1]), p=np.nan))
Rr = pd.DataFrame(res)
def meta(d):
    z = np.arctanh(d.rho.clip(-.999, .999)); w = d.n - 3; zf = (w * z).sum() / w.sum(); Q = (w * (z - zf) ** 2).sum()
    tau2 = max(0, (Q - (len(d) - 1)) / (w.sum() - (w ** 2).sum() / w.sum())); wr = 1 / (1 / w + tau2); zr = (wr * z).sum() / wr.sum(); se = wr.sum() ** -.5
    return dict(k=len(d), rho_random=np.tanh(zr), ci_lo=np.tanh(zr - 1.96 * se), ci_hi=np.tanh(zr + 1.96 * se), p=2 * (1 - norm.cdf(abs(zr / se))), I2=max(0, (Q - (len(d) - 1)) / Q) if Q > 0 else 0,
                frac_pos=(d.rho > 0).mean())
MT = pd.DataFrame({m: meta(g) for m, g in Rr.groupby("metric")}).T.round(3); MT.to_csv(f"{R}/pchembl_meta_analysis.csv"); Rr.to_csv(f"{R}/pchembl_per_uniprot_correlations.csv", index=False)
say("random-effects meta-analysis of per-UniProt Spearman correlations with pChEMBL:\n" + MT.to_string())

# ---------------- plots
fig, ax = plt.subplots(2, 2, figsize=(12, 9))
ax[0, 0].hist([C.water_B[C.interaction == k].clip(upper=100) for k in ["direct_water_metal", "direct_water_protein", "direct_water_only", "proximal_no_polar_contact"]], bins=50, density=True, stacked=False, histtype="step",
              label=["via metal", "via protein", "ligand only", "proximal"]); ax[0, 0].set(xlabel="water B (A$^2$), clipped at 100", title="Ligand-site water B by interaction"); ax[0, 0].legend()
ax[0, 1].hexbin(C.resolution, C.B_ratio.clip(upper=5), gridsize=40, bins="log", cmap="viridis"); ax[0, 1].set(xlabel="resolution (A)", ylabel="water B / mean protein B", title="Normalised water B vs resolution")
ax[1, 0].hist(S.frac_present, bins=40); ax[1, 0].set(xlabel="fraction of pocket structures with a water at the site", ylabel="sites", title="Water-site occupancy")
for i, (m, g) in enumerate(Rr.groupby("metric")): ax[1, 1].scatter(g.rho, np.full(len(g), i) + np.random.uniform(-.2, .2, len(g)), s=6, alpha=.5)
ax[1, 1].set(yticks=range(Rr.metric.nunique()), yticklabels=sorted(Rr.metric.unique()), xlabel="Spearman rho with pChEMBL (per UniProt)", title="pChEMBL vs water metrics"); ax[1, 1].axvline(0, c="k", lw=.5)
plt.tight_layout(); plt.savefig(f"{R}/summary_plots.png", dpi=150)
open(f"{R}/summary.txt", "w", encoding="utf-8").write("\n".join(out))
