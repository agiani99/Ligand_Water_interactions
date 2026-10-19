"""Pilot step 2: download mmCIF, find waters within 5 A of each ligand instance, classify interaction, record B-factors."""
import os, sys, requests, gemmi, numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor

DIRECT, PROX, METAL = 3.5, 5.0, 2.8
BB = {"N", "CA", "C", "O", "OXT"}
METALS = {"MG","ZN","MN","CA","FE","FE2","CO","NI","CU","CD","NA","K","HG","CU1","CU2","MO","W"}
os.makedirs("cif", exist_ok=True)
df = pd.read_csv("entries_ligands.csv")

def get(pdb):
    f = f"cif/{pdb}.cif"
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        r = requests.get(f"https://files.rcsb.org/download/{pdb}.cif", timeout=60); r.raise_for_status()
        open(f, "wb").write(r.content)
    return f

def analyse(pdb, ligs):
    st = gemmi.read_structure(get(pdb)); st.setup_entities(); st.remove_hydrogens()
    model = st[0]; ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
    prot = [a.b_iso for ch in model for r in ch if gemmi.find_tabulated_residue(r.name) and gemmi.find_tabulated_residue(r.name).is_amino_acid() for a in r]
    bprot = float(np.mean(prot)) if prot else np.nan
    out = []
    for ch in model:
        for res in ch:
            if res.name not in ligs: continue
            lig_pos = [(a, a.pos) for a in res]
            waters = {}
            for a, p in lig_pos:
                for m in ns.find_atoms(p, '\0', radius=PROX):
                    cra = m.to_cra(model)
                    if cra.residue.name != "HOH": continue
                    im = st.cell.find_nearest_image(p, m.pos, gemmi.Asu.Any); d = im.dist()
                    key = (cra.chain.name, str(cra.residue.seqid), cra.atom.altloc)
                    w = waters.setdefault(key, dict(cra=cra, pos=st.cell.find_nearest_pbc_position(p, cra.atom.pos, m.image_idx), dmin=99, dpol=99, sym=im.sym_idx != 0))
                    w["dmin"] = min(w["dmin"], d)
                    if a.element.name in ("N", "O"): w["dpol"] = min(w["dpol"], d)
            for key, w in waters.items():
                if w["dmin"] > PROX: continue
                wp = w["pos"]; pd_, pname, pbb, md, mname = 99, "", None, 99, ""
                for m in ns.find_atoms(wp, '\0', radius=4.0):
                    cra = m.to_cra(model); rn = cra.residue.name; e = cra.atom.element.name
                    if rn == "HOH" or cra.residue is res: continue
                    d = st.cell.find_nearest_image(wp, m.pos, gemmi.Asu.Any).dist()
                    ri = gemmi.find_tabulated_residue(rn)
                    if ri and ri.is_amino_acid() and e in ("N", "O") and d < pd_:
                        pd_, pname, pbb = d, f"{rn}{cra.residue.seqid.num}:{cra.atom.name}", cra.atom.name in BB
                    elif cra.atom.element.is_metal or rn in METALS:
                        if d < md: md, mname = d, f"{rn}{cra.residue.seqid.num}"
                cra = w["cra"]
                direct = w["dpol"] <= DIRECT
                if direct and md <= METAL: itype = "direct_water_metal"
                elif direct and pd_ <= DIRECT: itype = "direct_water_protein"
                elif direct: itype = "direct_water_only"
                else: itype = "proximal_no_polar_contact"
                out.append(dict(pdb=pdb, ligand=res.name, lig_chain=ch.name, lig_res=f"{res.name}{res.seqid.num}", water=f"{cra.chain.name}{cra.residue.seqid.num}",
                    water_occ=cra.atom.occ, water_B=cra.atom.b_iso, protein_meanB=bprot, B_ratio=cra.atom.b_iso/bprot if bprot else np.nan,
                    d_lig_min=round(w["dmin"],2), d_lig_polar=round(w["dpol"],2) if w["dpol"]<99 else np.nan, interaction=itype,
                    protein_partner=pname, d_protein=round(pd_,2) if pd_<99 else np.nan, partner_is_backbone=pbb,
                    metal=mname, d_metal=round(md,2) if md<99 else np.nan, symmetry_water=w["sym"],
                    wx=wp.x, wy=wp.y, wz=wp.z))
    return out

def job(args):
    pdb, ligs = args
    try: return analyse(pdb, set(ligs))
    except Exception as e: print("FAIL", pdb, e, file=sys.stderr); return []

groups = [(p, g.ligand.unique().tolist()) for p, g in df.groupby("pdb")]
with ThreadPoolExecutor(8) as ex: res = list(ex.map(job, groups))
w = pd.DataFrame([r for rs in res for r in rs])
w = w.merge(df.drop(columns=["seq","smiles"]), on=["pdb","ligand"], how="left")
w.to_csv("water_contacts.csv", index=False)
print(len(w), "ligand-water pairs;", w.pdb.nunique(), "structures")
print(w.groupby(["target","interaction"]).size().unstack(fill_value=0))
