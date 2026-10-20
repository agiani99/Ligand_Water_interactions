"""Full run, step 1: all human X-ray PDB entries (res < 2.5 A) with ligands 250 < MW < 1000.

Writes (in ./data):
  entries_ligands.csv   one row per (entry, ligand) that passes the filters
  excluded_ligands.csv  rows removed by the exclusion rules (for audit)
  entity_map.csv        one row per polymer entity (entry, entity_id, UniProt, organism)
  target_counts.csv     UniProt accessions ranked by number of entries with a kept ligand
Nothing is downloaded except API responses; no structure files are written.
"""
import os, time, requests, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

MW_MIN, MW_MAX, RES_MAX = 250, 1000, 2.5
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "data"); os.makedirs(OUT, exist_ok=True)
EXCL = {l.split("#")[0].strip().upper() for l in open(os.path.join(HERE, "exclude_ligands.txt")) if l.split("#")[0].strip()}

def post(url, payload, tries=5):
    for i in range(tries):
        try:
            r = requests.post(url, json=payload, timeout=180); r.raise_for_status(); return r.json()
        except Exception as e:
            err = e; time.sleep(3 * (i + 1))
    raise err

def search():
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entity_source_organism.taxonomy_lineage.name", "operator": "exact_match", "value": "Homo sapiens"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "exptl.method", "operator": "exact_match", "value": "X-RAY DIFFRACTION"}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less", "value": RES_MAX}},
        {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.nonpolymer_entity_count", "operator": "greater", "value": 0}}]},
        "return_type": "entry", "request_options": {"return_all_hits": True}}
    return [x["identifier"] for x in post("https://search.rcsb.org/rcsbsearch/v2/query", q).get("result_set", [])]

GQL = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
 rcsb_entry_info{resolution_combined}
 polymer_entities{rcsb_polymer_entity_container_identifiers{entity_id reference_sequence_identifiers{database_accession database_name}}
   rcsb_entity_source_organism{ncbi_scientific_name} entity_poly{rcsb_entity_polymer_type}}
 nonpolymer_entities{rcsb_nonpolymer_entity_container_identifiers{nonpolymer_comp_id}
   nonpolymer_comp{chem_comp{formula_weight type} rcsb_chem_comp_descriptor{SMILES InChIKey}}}}}"""

_hac = {}
def heavy_atoms(cid, smi):
    if cid not in _hac:
        m = Chem.MolFromSmiles(smi, sanitize=False) if smi else None
        _hac[cid] = m.GetNumHeavyAtoms() if m is not None else None
    return _hac[cid]

if __name__ == "__main__":
    ids = search(); print(len(ids), "entries match (human, X-ray, <2.5 A, >=1 non-polymer entity)", flush=True)
    rows, excl, ents = [], [], []
    for i in range(0, len(ids), 100):
        for e in post("https://data.rcsb.org/graphql", {"query": GQL, "variables": {"ids": ids[i:i+100]}})["data"]["entries"] or []:
            pdb, res = e["rcsb_id"], (e["rcsb_entry_info"]["resolution_combined"] or [None])[0]
            for pe in e["polymer_entities"] or []:
                c = pe["rcsb_polymer_entity_container_identifiers"]
                ups = [s["database_accession"] for s in (c.get("reference_sequence_identifiers") or []) if s["database_name"] == "UniProt"]
                org = ";".join(sorted({o["ncbi_scientific_name"] for o in (pe.get("rcsb_entity_source_organism") or []) if o.get("ncbi_scientific_name")}))
                ents.append(dict(pdb=pdb, entity_id=c["entity_id"], uniprot=ups[0] if ups else "", organism=org, polymer_type=(pe.get("entity_poly") or {}).get("rcsb_entity_polymer_type")))
            for ne in e["nonpolymer_entities"] or []:
                cid = ne["rcsb_nonpolymer_entity_container_identifiers"]["nonpolymer_comp_id"]
                comp = ne["nonpolymer_comp"]; cc = comp["chem_comp"]; d = comp.get("rcsb_chem_comp_descriptor") or {}
                mw = cc["formula_weight"]
                if mw is None or not (MW_MIN < mw < MW_MAX): continue
                row = dict(pdb=pdb, resolution=res, ligand=cid, mw=mw, heavy_atoms=heavy_atoms(cid, d.get("SMILES")), comp_type=cc.get("type"), inchikey=d.get("InChIKey"), smiles=d.get("SMILES"))
                why = "listed_additive" if cid.upper() in EXCL else ("saccharide" if "SACCHARIDE" in (cc.get("type") or "").upper() else "")
                (excl if why else rows).append({**row, "excluded_reason": why} if why else row)
        if (i // 100) % 50 == 0: print(f"  graphql {min(i+100, len(ids))}/{len(ids)}", flush=True)
    L, X, E = pd.DataFrame(rows), pd.DataFrame(excl), pd.DataFrame(ents)
    L.to_csv(f"{OUT}/entries_ligands.csv", index=False); X.to_csv(f"{OUT}/excluded_ligands.csv", index=False); E.to_csv(f"{OUT}/entity_map.csv", index=False)
    print(f"\nkept: {len(L)} entry-ligand rows, {L.pdb.nunique()} entries, {L.ligand.nunique()} distinct ligands")
    print(f"excluded: {len(X)} rows ({X.excluded_reason.value_counts().to_dict() if len(X) else {}})")
    human = E[E.organism.str.contains("Homo sapiens", na=False) & (E.uniprot != "")]
    T = human[human.pdb.isin(L.pdb)].groupby("uniprot").pdb.nunique().sort_values(ascending=False).rename("n_entries").reset_index()
    T.to_csv(f"{OUT}/target_counts.csv", index=False)
    print(f"{len(T)} human UniProt targets; with >=5 entries: {(T.n_entries >= 5).sum()}; top 10:\n{T.head(10).to_string(index=False)}")
