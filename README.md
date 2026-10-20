# Water_X: ligand-bound waters in human PDB structures

Analysis of crystallographic waters near co-crystallised ligands across **all human X-ray structures in the PDB** at resolution < 2.5 Å. It asks four questions:
1. How are ligand-site water B-factors distributed, and does a generic "textbook" summary of them (the Gemini answer, `Gemini_ansqer.md`) hold up?
2. How do B-factors depend on resolution and on the kind of contact?
3. How often does a ligand displace a conserved water?
4. Does displacing a water, or interacting with one, relate to ligand potency (pChEMBL)?

**Short answer:** pocket waters are mostly well-ordered, with the lowest B for metal- and protein-bridged waters, and absolute B cutoffs depend on resolution. Ligands displace conserved waters often (6,153 events). But there is **no detectable potency advantage for water-displacing ligands** once ligand size is controlled, in either a correlation analysis or a matched-pair test.

## Dataset and definitions

| Filter | Value |
|---|---|
| Organism / method | Homo sapiens, X-ray |
| Resolution | < 2.5 Å |
| Ligand MW | 250 < MW < 1000 (RCSB formula weight; heavy atoms from RDKit) |
| Excluded | glycans (by chemical type) and listed additives/detergents (`full/exclude_ligands.txt`) |
| Entries matched / kept | 40,672 / 23,052 (26,126 entry-ligand rows, 15,085 distinct ligands) |

- **Direct:** water within 3.5 Å of a ligand N or O.
- **Via protein / via metal:** a direct water that is also within 3.5 Å of a protein N or O, or within 2.8 Å of a metal.
- **Proximal:** no polar contact, but within 5 Å of any ligand atom. The minimum distance is reported.
- **B ratio:** water B divided by the mean protein B of the same structure (normalises across structures).
- **Water site:** waters from all superimposed structures of one UniProt, clustered at 1.2 Å. Pockets are superimposed on the Cα atoms within 12 Å of the ligand.
- **Conserved site:** a site with a water in ≥ 50% of pocket structures (≥ 10 structures).
- **Displaced:** the site has no water and a ligand atom lies within 2.5 Å of it.
- **pChEMBL:** local ChEMBL 37, binding assays with confidence score 8 or 9, joined by InChIKey and matched on the UniProt of the target.

| Quantity | Value |
|---|---|
| Ligand-water pairs (≤ 5 Å, excluding symmetry waters) | 373,812 in 22,033 entries, 1,821 UniProts |
| Water sites (UniProts with ≥ 5 structures) | 17,482 in 644 UniProts; 3,081 conserved |
| Pockets superposed | 22,860 fitted (pocket RMSD ≤ 2.0 Å); 546 skipped |
| Ligand-water displacement events at conserved sites | 6,153 in 4,244 structures, 336 UniProts, 3,143 ligands |
| Ligands with same-target pChEMBL | 6,107 (UniProt, ligand) pairs in 648 UniProts; 109 UniProts with ≥ 15 ligands |

## Findings vs the Gemini textbook answer

| Gemini claim | Result | Verdict |
|---|---|---|
| Ligand-site waters are "well-ordered", ~20–35 Å² | Median 32.4 Å² (IQR 24–42). 41% of pocket waters fall in 20–35 Å², but 39% are in 35–60 Å² | **Partly supported** |
| Anchored, tightly H-bonded waters have the lowest B | Median B: metal-bridged 23 Å², protein-bridged 29, proximal 33, ligand-only 36. Waters H-bonded to the ligand have a lower resolution-corrected B ratio at 61% of 924 sites (median difference −0.043, p = 2e-17), a real but small effect | **Supported** |
| Waters < 20 Å² are structural "linchpins" | 16% of pocket waters are < 20 Å², far too many to all be linchpins. Low B alone doesn't identify them | **Overstated** |
| > 45–60 Å² is mobile or bulk-like | The 95th percentile is 57 Å² and only 3% of pocket waters exceed 60 Å². Pocket waters are rarely bulk-like | **Weak fit** |
| > 80 Å² is a "red flag" | 0.3% of pocket waters exceed 80 Å² | **Real but rare** |
| Fixed absolute cutoffs (20/35/60 Å²) | Median B rises from 18 Å² (≤ 1.0 Å) to 27, 32 and 38 Å² (2.0–2.5 Å) (Spearman ρ = 0.29). Median B per UniProt spans IQR 29–42 Å². The normalised B ratio falls with resolution value (ρ = −0.41) because poorly resolved structures keep only their best-ordered waters | **Not supported: use B ratio and resolution matching** |
| Displacing a bound water "dramatically" boosts affinity | No detectable effect once ligand size is controlled (see the two sections below) | **Not supported in this data** |

Overall: the broad ranking (metal- and protein-bridged waters are the most ordered) is supported. The fixed thresholds and the affinity claim are not.

## Water sites, resolution and displacement

- The fraction of pocket sites holding a water falls at poorer resolution (ρ = −0.16), and so does the number of displacement calls (ρ = −0.15).
- At conserved sites, 5.5% of pocket observations are called displaced.
- Top targets by number of displacement events: P01116 (KRAS), P56817 (BACE1), O60885 (BRD4), P41182 (BCL6), P15090 and P00918 (carbonic anhydrase II).

## pChEMBL vs water interaction

Random-effects meta-analysis of per-UniProt Spearman correlations (109 UniProts with ≥ 15 ligands):

| Metric | ρ | 95% CI | After controlling ligand size |
|---|---|---|---|
| Direct water contacts | 0.089 | 0.039–0.139 | 0.039 (−0.009 to 0.088) |
| Bridging waters | 0.069 | 0.020–0.117 | 0.028 (−0.020 to 0.075) |
| Displaced conserved waters | 0.038 | −0.012 to 0.088 | −0.030 (−0.079 to 0.019) |
| B ratio of H-bonded waters | −0.022 | −0.066 to 0.021 | −0.030 (−0.075 to 0.014) |

Raw correlations for contact counts are small and positive, but they vanish once ligand heavy-atom count is controlled. Within a target, pChEMBL tracks ligand size strongly (ρ = 0.37). Pooling targets without this adjustment would be misleading because of target identity and ligand size.

## Matched-pair test: displaced vs retained water

Correlation can't fully separate water effects from ligand size, so I compared pairs of **different ligands of the same target at the same conserved water site**. Ligand A displaces the water, ligand B keeps it, the two are structurally similar (Morgan-2 Tanimoto), were crystallised at similar resolution (|Δresolution| ≤ 0.5 Å), and both have a pChEMBL for that target. The effect is pChEMBL(A) − pChEMBL(B), so a positive value means the water-displacing ligand is more potent. A ligand pair counts once per target, and confidence intervals use a bootstrap over targets, since pairs within one target are not independent.

Main setting (Tanimoto ≥ 0.6, |Δres| ≤ 0.5 Å, sites occupied in ≥ 50% of pocket structures):

| | Value |
|---|---|
| Matched ligand pairs | 280 in 66 targets |
| Mean ΔpChEMBL (displacer − retainer) | +0.05 (95% CI −0.10 to +0.19) |
| Adjusted for heavy-atom and resolution difference | −0.05 (−0.18 to +0.07) |
| Fraction of pairs where the displacer is more potent | 54% |
| Mean heavy-atom difference | +1.1 (displacers are larger) |

Sensitivity (full table in `full/data/results/matched_pairs_sensitivity.csv`):
- **Looser similarity (Tanimoto ≥ 0.5, 599 pairs, 93 targets):** the raw difference is +0.13 (CI +0.01 to +0.26). After adjusting for heavy-atom difference it is −0.03 (−0.12 to +0.10), so the apparent gain is explained by ligand size.
- **Stricter similarity (≥ 0.7 or ≥ 0.8) and tighter resolution (≤ 0.3 Å):** the estimates stay between −0.06 and +0.09, all with CIs spanning zero.
- **Same-size pairs only (|Δ heavy atoms| ≤ 3 or ≤ 1):** −0.01 (CI −0.15 to +0.11, 210 pairs) and −0.01 (−0.16 to +0.13, 140 pairs).

Interpretation: there is **no detectable potency advantage for ligands that displace a conserved crystallographic water**. The data exclude effects larger than about 0.2–0.3 log units (roughly 2-fold) at this sample size. They cannot rule out a small effect, and they do not test water displacement in specific series where the water is known to be energetically unfavourable.

## Caveats

- Ligand binding is assigned to the protein chain with the most contacts. Only the first ligand copy per entry is used in the site analysis.
- Displacement is geometric (a site with no water and a ligand atom within 2.5 Å) and cannot separate a true displacement from an unmodelled water.
- Water B is normalised by the structure-wide protein mean B.
- Cutoffs (3.5 / 2.8 / 5.0 Å) are conventional defaults, not tuned.
- pChEMBL values come from different assays and papers and exist for only about 6,100 (UniProt, ligand) pairs out of 15,000 ligands.
- Ligand similarity in the matched-pair test is whole-molecule fingerprint similarity, not a matched molecular pair on one substituent.

## Reproducing

Python 3.12 with gemmi, rdkit, pandas, numpy, scipy, matplotlib, requests, biopython. Run the scripts in `full/` in order from that folder. Steps 1–2 need internet access: structures are streamed from RCSB in memory and never written to disk. Step 4 reads a local ChEMBL SQLite file (path set at the top of `04_chembl.py`).

| Script | Output (in `full/data/`) |
|---|---|
| `01_select.py` | `entries_ligands.csv`, `excluded_ligands.csv`, `entity_map.csv`, `target_counts.csv` |
| `02_contacts.py` | `contacts/` (per-chunk ligand-water contacts and pocket data; resumable) |
| `03_sites.py` | `sites/` (fits, water sites, site-by-structure table) |
| `04_chembl.py` | `chembl/ligand_pchembl.csv` |
| `05_stats.py` | `results/` (`summary.txt`, `summary_plots.png`, correlation tables) |
| `06_matched_pairs.py` | `results/matched_pairs_*.csv`, `matched_pairs.png` |
