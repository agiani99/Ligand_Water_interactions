# Water_X: ligand-bound waters in human PDB structures (pilot)

Pilot analysis of crystallographic waters near co-crystallised ligands, run on three human targets:
**CDK2** (P24941), **HSP90α** (P07900) and **carbonic anhydrase II** (P00918).
The pilot also tests how well a generic "textbook" summary of water B-factors (the Gemini answer, `Gemini_ansqer.md`) holds up against real data.

## Dataset

| Filter | Value |
|---|---|
| Organism / method | Homo sapiens, X-ray |
| Resolution | < 2.5 Å |
| Ligand MW | 250 < MW < 1000 (RCSB formula weight; heavy atoms from RDKit) |
| Structures | 1,294 (CDK2 325, HSP90α 248, CA2 721) |
| Ligand-water pairs (water ≤ 5 Å from a ligand atom) | 14,265 |

Definitions used:
- **Direct:** water within 3.5 Å of a ligand N or O.
- **Via protein / via metal:** a direct water that is also within 3.5 Å of a protein N or O, or within 2.8 Å of a metal.
- **Proximal:** no polar contact, but within 5 Å of any ligand atom. The minimum distance is reported.
- **B ratio:** water B divided by the mean protein B of the same structure (normalises across structures).
- **Water site:** waters from all superimposed structures of a target, clustered at 1.2 Å.
- **Conserved site:** a site with a water in ≥ 50% of pocket structures.
- **Displaced:** the site has no water and a ligand atom lies within 2.5 Å of it.

## Findings vs the Gemini textbook answer

| Gemini claim | Pilot result | Verdict |
|---|---|---|
| Ligand-site waters are "well-ordered", ~20–35 Å² | Median 28.6 Å² (IQR 21–37). 49% of pocket waters fall in 20–35 Å² | **Supported** |
| Anchored, tightly H-bonded waters have the lowest B | Lowest median B for metal-bridging waters (19 Å², n = 105), then protein-bridging (25 Å², n = 3,148). Proximal-only waters 29 Å², ligand-only 33 Å² | **Supported** |
| Waters < 20 Å² are structural "linchpins" | 21% of pocket waters are < 20 Å², too many to all be linchpins. Low B alone doesn't identify them | **Overstated** |
| > 45–60 Å² is mobile or bulk-like | Only ~3% of pocket waters are > 45 Å² (95th percentile 52 Å²). Pocket waters are rarely bulk-like | **Weak fit** |
| > 80 Å² is a "red flag" | 0.3% of pocket waters exceed 80 Å² (max 154 Å²) | **Real but rare** |
| Fixed absolute cutoffs (20/35/60 Å²) | Median B shifts with resolution (27 Å² at ≤ 1.5 Å, 29 Å² at 1.5–2.0 Å, 36 Å² at 2.0–2.5 Å) and by target (CA2 27, HSP90α 27, CDK2 36 Å²). Cutoffs would mislabel whole structures | **Not supported: use B ratio and resolution matching** |
| Displacing a bound water "dramatically" boosts affinity | No consistent pChEMBL effect (see below) | **Not supported in this data** |

Overall: the broad ranking (metal- and protein-bridged waters are the most ordered) is supported. The fixed thresholds and the affinity claim are not.

## Other results

**Resolution**
- Water B rises with resolution value (Spearman ρ = 0.18, n = 11,025 waters at conserved and non-conserved sites).
- B ratio falls (ρ = −0.47).
- The fraction of conserved sites with a modelled water falls at poorer resolution (ρ = −0.30).
- Poorly resolved structures keep only their best-ordered waters, so the effect depends on how B is normalised.

**Conserved sites and displacement**
- 564 sites were found, 53 of them conserved.
- 371 ligand-water displacement events at conserved sites, in 311 structures.
- Per target: CA2 180 events (155 structures), CDK2 85 (82), HSP90α 106 (74).

**pChEMBL vs water (step 10)**
- Source: ChEMBL binding assays with confidence score 8 or 9, joined to ligands by InChIKey.
- 603 of 1,003 ligands are in ChEMBL, and 366 have a pChEMBL value for the same target (CA2 187, CDK2 90, HSP90α 88).
- No consistent correlation was found, either per target or after correcting for ligand heavy-atom count:

| Metric | CA2 ρ | CDK2 ρ | HSP90α ρ |
|---|---|---|---|
| Direct water contacts | −0.02 | −0.05 | 0.08 |
| Bridging waters | 0.02 | 0.04 | 0.07 |
| Displaced conserved waters | **−0.18** (p = 0.014) | −0.04 | 0.11 |

- The CA2 trend is weak. Ligands that displace a conserved water have a median pChEMBL of 7.30 vs 7.77. The CA2 result does not replicate in CDK2 or HSP90α.
- A pooled correlation between mean water B ratio and pChEMBL (ρ = 0.22, all targets) is driven by target identity. Do not interpret it.

## Caveats

- Only three targets, so the numbers are not general. The PDB is biased toward a few targets.
- Displacement is inferred geometrically, so an unmodelled water can't be told apart from a true displacement.
- Only the first copy of each ligand per structure is used in the site analysis.
- About 15 CDK2 structures in a different conformation, and 2 HSP90α structures, failed whole-chain superposition and were skipped.
- Cutoffs (3.5 / 2.8 / 5.0 Å) are conventional defaults, not tuned.
- pChEMBL coverage is limited to ligands present in ChEMBL (about 36% of PDB ligands here).
- Gemini's literature citations were not verified.

## Reproducing

Python 3.12 with gemmi, rdkit, pandas, numpy, scipy, matplotlib, requests. Run the scripts in `pilot/` in order from that folder. They need internet access and download about 1,300 mmCIF files into `pilot/cif/`.

| Script | Output |
|---|---|
| `01_fetch_entries.py` | `entries_ligands.csv`: structures and ligands (RCSB) |
| `02_water_contacts.py` | `water_contacts.csv`: per ligand-water pair: interaction, distances, B |
| `03_water_sites.py` | `obs.pkl`: superimposed pocket waters and ligands |
| `04_analysis.py` | `water_sites.csv`, `site_structure_table.csv`, `water_B_direct_vs_nocontact.csv`, `water_B_plots.png` |
| `05_chembl.py` | `ligand_pchembl.csv` (ChEMBL API) |
| `06_pchembl_vs_water.py` | `pdb_ligand_pchembl_water.csv`, `ligand_level_pchembl_water.csv` |

## Full-scale run (`full/`)

Same conditions as the pilot, extended to **all** human X-ray structures: 40,672 matched, 23,052 kept with a ligand 250 < MW < 1000 after excluding glycans and listed additives (`full/exclude_ligands.txt`). Structures are streamed from RCSB in memory and never written to disk (`02_contacts.py`).

| Quantity | Value |
|---|---|
| Ligand-water pairs (≤ 5 Å) | 373,812 (non-symmetry waters) in 22,033 entries, 1,821 UniProts |
| Water sites (UniProt groups with ≥ 5 structures) | 17,482 in 644 UniProts; 3,081 conserved (≥ 50% occupancy, ≥ 10 structures) |
| Ligand-water displacement events at conserved sites | 6,153 in 4,244 structures, 336 UniProts, 3,143 ligands |
| Pockets superposed | 22,860 fitted (pocket RMSD ≤ 2.0 Å); 546 skipped |
| Ligands with same-target pChEMBL (local ChEMBL 37, binding assays, confidence 8/9) | 6,107 (UniProt, ligand) pairs in 648 UniProts; 109 UniProts with ≥ 15 ligands |

**Water B-factors vs the Gemini bands.** Median 32.4 Å² (IQR 24–42). The fractions are below 20 Å² 16%, 20–35 Å² 41%, 35–60 Å² 39%, 60–80 Å² 3%, above 80 Å² 0.3%. Most pocket waters therefore sit in the "well-ordered to mobile" range and few are extreme, but the 20–35 band is not special. Median B by interaction: metal-bridged 23, protein-bridged 29, proximal 33, ligand-only 36 Å². The ordering supports the "anchored means lower B" idea. By resolution, median B goes from 18 Å² (≤ 1.0 Å) to 38 Å² (2.0–2.5 Å) (Spearman ρ = 0.29), so absolute cutoffs are resolution-dependent. Normalised B ratio falls with resolution value (ρ = −0.41).

**Direct contact vs none.** Within 924 sites with enough data, waters H-bonded to the ligand have a lower resolution-corrected B ratio at 61% of sites (median difference −0.043, Wilcoxon p = 2e-17). The effect is real but small.

**Resolution and displacement.** The fraction of pocket sites holding a water falls at poorer resolution (ρ = −0.16), and so does the number of displacement calls (ρ = −0.15). At conserved sites 5.5% of pocket observations are called displaced. The top targets by events are P01116 (KRAS), P56817 (BACE1), O60885 (BRD4), P41182 (BCL6), P15090 and P00918 (CA2).

**pChEMBL vs water interaction (step 10).** Random-effects meta-analysis of per-UniProt Spearman correlations (109 UniProts):

| Metric | ρ | 95% CI | After controlling ligand size |
|---|---|---|---|
| Direct water contacts | 0.089 | 0.039–0.139 | 0.039 (−0.009 to 0.088) |
| Bridging waters | 0.069 | 0.020–0.117 | 0.028 (−0.020 to 0.075) |
| Displaced conserved waters | 0.038 | −0.012 to 0.088 | −0.030 (−0.079 to 0.019) |
| B ratio of H-bonded waters | −0.022 | −0.066 to 0.021 | −0.030 (−0.075 to 0.014) |

Raw correlations for contact counts are small and positive, but they vanish once ligand heavy-atom count is controlled. Within a target, pChEMBL tracks ligand size strongly (ρ = 0.37). Conclusion: **no detectable relationship between water interactions or displacement and potency in this dataset, beyond ligand size.** The pilot's weak CA2 trend does not hold at scale.

Full-scale caveats: ligand binding is assigned to the protein chain with the most contacts; only the first ligand copy per entry is used for sites; displacement is geometric (a site with no water and a ligand atom within 2.5 Å) and cannot separate true displacement from an unmodelled water; the structure-wide protein mean B is used for normalisation; pChEMBL coverage is about 6,100 pairs for 15,000 ligands.

Run order in `full/`: `01_select.py` → `02_contacts.py` → `03_sites.py` → `04_chembl.py` → `05_stats.py`. Outputs: `full/data/results/` (`summary.txt`, `summary_plots.png`, CSV tables).

## Next steps

- Validate displacement calls against known cases (HSP90, CDK2, BRD4, BACE1, KRAS).
- Stratify results by protein family and ligand chemotype.
- Compare same-target pairs with and without a displaced water at matched resolution and matched ligand series, which is stronger than correlation.
