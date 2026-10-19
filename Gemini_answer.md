In macromolecular X-ray crystallography, the B-factor (or temperature factor) measures how much an atom vibrates or fluctuates around its average position. For water molecules trapped inside or coating a protein structure, their average B-factor values are incredibly telling, signaling whether a water molecule is a permanent structural fixture or just a passing, transient solvent molecule. [1] (https://pmc.ncbi.nlm.nih.gov/articles/PMC6899928/), [2] (https://www.youtube.com/watch?v=pRdmXs_KbE4&t=155), [3] (https://www.proteinstructures.com/pdb-file-format/), [4] (https://pubs.acs.org/jcisd8/article/55/4/771/941552/Evidence-of-Water-Molecules-amp-xe5f8-A), [5] (https://pubs.acs.org/jcisd8/article/47/2/668/3477884/Analysis-of-Ligand-Bound-Water-Molecules-in-High)
The average B-factor values that reveal interesting behavior in water molecules follow a distinct structural hierarchy:

1. Highly Conserved "Structural" Water: < 20 Å²

When a water molecule's average B-factor drops below 20 Å², it is behaving more like a rigid stone than a liquid fluid. [1] (https://www.youtube.com/watch?v=pRdmXs_KbE4&t=155)
• Why it's interesting: These values are frequently lower than the average B-factor of the surrounding protein atoms. This occurs when a water molecule is fully buried within the protein core and tightly locked down by 3 or 4 hydrogen bonds. [1] (https://pmc.ncbi.nlm.nih.gov/articles/PMC6899928/)
• Significance: These are structural "linchpins" that bridge separate pieces of the protein backbone together. They are so stable that they are highly conserved across evolutionarily related proteins. [1] (https://pmc.ncbi.nlm.nih.gov/articles/PMC6899928/)

2. High-Affinity "Ligand-Bound" Water: ~20 to 35 Å²

An average B-factor in this moderate, controlled range indicates a water molecule that is well-ordered and anchored to a specific binding site, pocket, or active site. [1] (https://pubs.acs.org/jcisd8/article/47/2/668/3477884/Analysis-of-Ligand-Bound-Water-Molecules-in-High), [2] (https://www.youtube.com/watch?v=pRdmXs_KbE4&t=155)
• Why it's interesting: These waters sit right at the functional interface. They are highly relevant to drug design.
• Significance: If a pharmaceutical drug candidate can be designed to mimic the exact geometry of this water molecule, it can displace the water, releasing it into bulk solvent. This release increases entropy, dramatically boosting the drug's binding affinity. [1] (https://pubs.rsc.org/cp/article/21/33/18149/639174/Are-crystallographic-B-factors-suitable-for), [2] (https://pubs.acs.org/jcisd8/article/47/2/668/3477884/Analysis-of-Ligand-Bound-Water-Molecules-in-High), [3] (https://pubs.acs.org/chreay/article/119/3/1626/1372803/Utility-of-B-Factors-in-Protein-Science), [4] (https://pubs.acs.org/jacsat/article/129/9/2577/3525810/Classification-of-Water-Molecules-in-Protein)

3. Mobile "Bulk-Like" Water: > 45 Å² to 60 Å²

When water molecules display an average B-factor above 45 Å², they are exhibiting rapid movement, structural flexibility, or partial occupancy (meaning the water is only present in some of the protein crystals in the lattice). [1] (https://chemistry.stackexchange.com/questions/156289/how-to-identify-water-mistakes-in-protein-crystal-structures-are-b-factors-an), [2] (https://onlinelibrary.wiley.com/doi/10.1002/jcc.23170), [3] (https://www.youtube.com/watch?v=pRdmXs_KbE4&t=155)
• Why it's interesting: This marks the boundary where water transitions from behaving like a localized part of the macromolecule to behaving like free, chaotic bulk water. [1] (https://pubmed.ncbi.nlm.nih.gov/29471780/), [2] (https://pubs.acs.org/jcisd8/article/47/2/668/3477884/Analysis-of-Ligand-Bound-Water-Molecules-in-High)
• Significance: When a water molecule forms fewer than 2 hydrogen bonds with a protein, it remains highly mobile and acts as a fluid cushion, facilitating loop movements and structural transitions necessary for the protein to change shape during a reaction. [1] (https://pubs.acs.org/chreay/article/119/3/1626/1372803/Utility-of-B-Factors-in-Protein-Science), [2] (https://pubs.acs.org/jcisd8/article/47/2/668/3477884/Analysis-of-Ligand-Bound-Water-Molecules-in-High)

4. The "Red Flag" Limit: > 80 Å² to 100 Å²

In structural biology, if a water molecule’s B-factor inflates above 80–100 Å², it crosses from "biologically interesting" to "crystallographically suspicious." [1] (https://proteopedia.org/index.php/Water_in_macromolecular_models)
• Why it's interesting: A B-factor this high means the electron density is so blurred out that the water molecule is practically invisible.
• Significance: It usually indicates that the crystal structure was over-refined or misassigned. Software or crystallographers often "guess" a water molecule into an unassigned blip of noise just to lower their statistical error (R-factor), even though an actual, static water molecule does not exist there. [1] (https://pmc.ncbi.nlm.nih.gov/articles/PMC5824579/), [2] (https://www.youtube.com/watch?v=pRdmXs_KbE4&t=155), [3] (https://pubmed.ncbi.nlm.nih.gov/29471780/), [4] (https://proteopedia.org/index.php/Water_in_macromolecular_models)

Directly Comparing Water Dynamics by B-Factor


Average B-Factor Range	Dynamic Behavior	Biological / Crystallographic Context
< 20 Å²	Immobile / Anchored	Buried deep in the core; bridges secondary structures; mimics rigid protein atoms.
20 – 35 Å²	Well-ordered / Functional	Active sites and binding pockets; top targets for displacement in drug design.
35 – 60 Å²	Mobile / Transient	Surface/hydration shell water; provides a lubricating fluid layer for protein movement.
> 80 Å²	Blurry / Disordered	Bulk solvent boundary or refinement errors; position is mostly speculative.

Are you looking at a specific PDB protein structure file, or are you exploring water B-factors for computational drug design or molecular dynamics simulations? Knowing your exact use case can help tailor the math or structural analysis guidelines.
