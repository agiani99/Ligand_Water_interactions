"""Pilot step 3: superimpose same-target structures, cluster pocket waters into sites,
and for each structure/site classify: present | displaced_by_ligand | absent."""
import os, sys, pickle, gemmi, numpy as np, pandas as pd
from collections import Counter
from scipy.spatial import cKDTree
from concurrent.futures import ThreadPoolExecutor

SITE_R, DISP_R, SHELL = 1.2, 2.5, 8.0
df = pd.read_csv("entries_ligands.csv")
cont = pd.read_csv("water_contacts.csv")

def observe(pdb, ligs, ref_poly=None):
    """first ligand instance per (pdb,ligand): chain used for fit, ligand atoms, waters<=SHELL. Returns list of dicts + fitted polymer."""
    st = gemmi.read_structure(f"cif/{pdb}.cif"); st.setup_entities(); st.remove_hydrogens()
    model = st[0]; ns = gemmi.NeighborSearch(model, st.cell, SHELL + 1).populate()
    obs = []; seen = set()
    for ch in model:
        for res in ch:
            if res.name not in ligs or res.name in seen: continue
            seen.add(res.name)
            cnt = Counter(); wat = {}
            for a in res:
                for m in ns.find_atoms(a.pos, '\0', radius=SHELL):
                    cra = m.to_cra(model)
                    if cra.residue is res or cra.residue.name == "HOH":
                        if cra.residue.name == "HOH":
                            im = st.cell.find_nearest_image(a.pos, m.pos, gemmi.Asu.Any)
                            if im.dist() <= SHELL:
                                k = (cra.chain.name, cra.residue.seqid.num, cra.atom.altloc)
                                wat.setdefault(k, (st.cell.find_nearest_pbc_position(a.pos, cra.atom.pos, m.image_idx), cra.atom.b_iso, cra.atom.occ))
                        continue
                    if m.image_idx == 0 and st.cell.find_nearest_image(a.pos, m.pos, gemmi.Asu.Any).dist() <= 5.0 and gemmi.find_tabulated_residue(cra.residue.name) \
                       and gemmi.find_tabulated_residue(cra.residue.name).is_amino_acid(): cnt[cra.chain.name] += 1
            if not cnt: continue
            fitch = model[cnt.most_common(1)[0][0]]; poly = fitch.get_polymer()
            if len(poly) < 50: continue
            obs.append(dict(pdb=pdb, ligand=res.name, lig=[(a.element.name, a.pos) for a in res], wat=list(wat.values()), poly=poly, model=model))
            break
    return obs

def run(target):
    d = df[df.target == target]
    ref_pdb = d.sort_values("resolution").iloc[0].pdb
    groups = [(p, set(g.ligand)) for p, g in d.groupby("pdb")]
    obs = []
    ref = observe(ref_pdb, set(d[d.pdb == ref_pdb].ligand))[0]; ref_poly = ref["poly"]; ptype = ref_poly.check_polymer_type()
    fixed = []
    for p, l in groups:
        try:
            for o in observe(p, l):
                sup = gemmi.calculate_superposition(ref_poly, o["poly"], ptype, gemmi.SupSelect.CaP)
                if sup.count < 50 or sup.rmsd > 3.0: print("skip fit", p, sup.count, round(sup.rmsd,2)); continue
                T = sup.transform
                o["lig"] = [(e, T.apply(q)) for e, q in o["lig"]]
                o["wat"] = [(T.apply(q), b, oc) for q, b, oc in o["wat"]]
                o["rmsd"] = sup.rmsd; del o["poly"], o["model"]; fixed.append(o)
        except Exception as e: print("FAIL", p, e, file=sys.stderr)
    return fixed

all_obs = {}
for t in df.target.unique():
    all_obs[t] = run(t); print(t, len(all_obs[t]), "structures fitted")
pickle.dump({t: [dict(pdb=o["pdb"], ligand=o["ligand"], rmsd=o["rmsd"],
    lig=[(e, (q.x,q.y,q.z)) for e,q in o["lig"]], wat=[((q.x,q.y,q.z),b,oc) for q,b,oc in o["wat"]]) for o in v] for t,v in all_obs.items()}, open("obs.pkl","wb"))
