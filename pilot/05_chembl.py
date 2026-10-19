"""Pilot step 5: PDB ligand -> ChEMBL compound (InChIKey) -> pChEMBL for same target (binding assays, confidence 8/9)."""
import requests, pandas as pd, time
B = "https://www.ebi.ac.uk/chembl/api/data"
TARGET = {"CDK2": "CHEMBL301", "HSP90AA1": "CHEMBL3880", "CA2": "CHEMBL205"}
ent = pd.read_csv("entries_ligands.csv")
comps = sorted(ent.ligand.unique())

# 1. InChIKeys from RCSB
GQL = "query($ids:[String!]!){chem_comps(comp_ids:$ids){chem_comp{id} rcsb_chem_comp_descriptor{InChIKey}}}"
key = {}
for i in range(0, len(comps), 100):
    r = requests.post("https://data.rcsb.org/graphql", json={"query": GQL, "variables": {"ids": comps[i:i+100]}}, timeout=120).json()
    for c in r["data"]["chem_comps"]: key[c["chem_comp"]["id"]] = (c["rcsb_chem_comp_descriptor"] or {}).get("InChIKey")
print(len(key), "ligand InChIKeys")

def get(path, **params):
    for _ in range(4):
        try:
            r = requests.get(f"{B}/{path}.json", params=params, timeout=90); r.raise_for_status(); return r.json()
        except Exception as e: err = e; time.sleep(3)
    raise err

# 2. InChIKey -> molecule_chembl_id
k2m = {}; keys = sorted({k for k in key.values() if k})
for i in range(0, len(keys), 40):
    j = get("molecule", molecule_structures__standard_inchi_key__in=",".join(keys[i:i+40]), limit=1000)
    for m in j["molecules"]: k2m[m["molecule_structures"]["standard_inchi_key"]] = m["molecule_chembl_id"]
print(len(k2m), "ligands found in ChEMBL")
mol = {c: k2m.get(k) for c, k in key.items()}

# 3. activities per target; assay filter (binding, confidence 8/9)
rows = []
for t, tid in TARGET.items():
    assays = set()
    off = 0
    while True:
        j = get("assay", target_chembl_id=tid, assay_type="B", confidence_score__in="8,9", limit=1000, offset=off)
        assays |= {a["assay_chembl_id"] for a in j["assays"]}
        if not j["page_meta"]["next"]: break
        off += 1000
    mids = sorted({m for c, m in mol.items() if m and c in set(ent[ent.target == t].ligand)})
    for i in range(0, len(mids), 40):
        off = 0
        while True:
            j = get("activity", target_chembl_id=tid, molecule_chembl_id__in=",".join(mids[i:i+40]), pchembl_value__isnull=False, limit=1000, offset=off)
            for a in j["activities"]:
                if a["assay_chembl_id"] in assays: rows.append(dict(target=t, molecule_chembl_id=a["molecule_chembl_id"], standard_type=a["standard_type"], pchembl=float(a["pchembl_value"])))
            if not j["page_meta"]["next"]: break
            off += 1000
    print(t, len(assays), "assays;", len(mids), "ligand compounds queried;", sum(r["target"] == t for r in rows), "activities")
A = pd.DataFrame(rows)
P = A.groupby(["target", "molecule_chembl_id"]).agg(pchembl=("pchembl", "median"), n_act=("pchembl", "size")).reset_index()
P["ligand"] = P.molecule_chembl_id.map({m: c for c, m in mol.items() if m})
P.to_csv("ligand_pchembl.csv", index=False); print(len(P), "target-ligand pChEMBL values")
