"""Pilot step 1: list PDB entries (human, res<2.5, X-ray) for pilot UniProts and their ligands (250<MW<1000)."""
import requests, json, pandas as pd
from rdkit import Chem

TARGETS = {"P24941": "CDK2", "P07900": "HSP90AA1", "P00918": "CA2"}
MW_MIN, MW_MAX, RES_MAX = 250, 1000, 2.5

def search(uniprot):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "operator": "exact_match", "value": uniprot}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entity_source_organism.taxonomy_lineage.name", "operator": "exact_match", "value": "Homo sapiens"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less", "value": RES_MAX}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.nonpolymer_entity_count", "operator": "greater", "value": 0}}]},
        "return_type": "entry", "request_options": {"return_all_hits": True}}
    r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=q, timeout=60)
    r.raise_for_status()
    return [x["identifier"] for x in r.json().get("result_set", [])]

GQL = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
 rcsb_entry_info{resolution_combined}
 polymer_entities{entity_poly{pdbx_seq_one_letter_code_can} rcsb_polymer_entity_container_identifiers{reference_sequence_identifiers{database_accession database_name}}}
 nonpolymer_entities{rcsb_nonpolymer_entity_container_identifiers{nonpolymer_comp_id} nonpolymer_comp{chem_comp{formula_weight} rcsb_chem_comp_descriptor{SMILES}}}}}"""

def fetch(ids):
    out = []
    for i in range(0, len(ids), 100):
        r = requests.post("https://data.rcsb.org/graphql", json={"query": GQL, "variables": {"ids": ids[i:i+100]}}, timeout=120)
        r.raise_for_status()
        out += r.json()["data"]["entries"]
    return out

rows = []
for up, name in TARGETS.items():
    ids = search(up)
    print(name, up, len(ids), "entries")
    for e in fetch(ids):
        res = e["rcsb_entry_info"]["resolution_combined"][0]
        seqs = [pe["entity_poly"]["pdbx_seq_one_letter_code_can"] for pe in e["polymer_entities"]
                if any(s["database_accession"] == up for s in (pe["rcsb_polymer_entity_container_identifiers"].get("reference_sequence_identifiers") or []))]
        seq = seqs[0] if seqs else ""
        for ne in e["nonpolymer_entities"] or []:
            cid = ne["rcsb_nonpolymer_entity_container_identifiers"]["nonpolymer_comp_id"]
            comp = ne["nonpolymer_comp"]
            mw = comp["chem_comp"]["formula_weight"]
            smi = (comp.get("rcsb_chem_comp_descriptor") or {}).get("SMILES")
            if mw is None or not (MW_MIN < mw < MW_MAX): continue
            hac = Chem.MolFromSmiles(smi).GetNumHeavyAtoms() if smi and Chem.MolFromSmiles(smi) else None
            rows.append(dict(target=name, uniprot=up, pdb=e["rcsb_id"], resolution=res, ligand=cid, mw=mw, heavy_atoms=hac, smiles=smi, seq=seq))
df = pd.DataFrame(rows)
df.to_csv("entries_ligands.csv", index=False)
print(df.groupby("target").agg(entries=("pdb", "nunique"), ligands=("ligand", "nunique"), rows=("pdb", "size")))
