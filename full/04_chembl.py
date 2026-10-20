"""Full run, step 4: pChEMBL for PDB ligands from the LOCAL ChEMBL SQLite (read-only).

PDB ligand -> InChIKey -> ChEMBL compound; target matched by UniProt accession of single-protein targets;
binding assays (assay_type 'B') with confidence_score 8 or 9; pChEMBL not null.
Output: data/chembl/ligand_pchembl.csv  (uniprot, ligand, pchembl median over measurements, n_act, n_assays)
"""
import os, sqlite3, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data"); os.makedirs(f"{DATA}/chembl", exist_ok=True)
DB = r"C:\Users\andrea.DESKTOP-26V6UN4\Documents\NEMO_AI\act_cliff\SQLITE\chembl_37_sqlite\chembl_37.db"

L = pd.read_csv(f"{DATA}/entries_ligands.csv").dropna(subset=["inchikey"]).drop_duplicates("ligand")[["ligand", "inchikey"]]
E = pd.read_csv(f"{DATA}/entity_map.csv", dtype=str).fillna("")
ups = sorted(set(E[(E.uniprot != "") & E.organism.str.contains("Homo sapiens")].uniprot))
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)           # read-only: the db has -wal/-shm files, never write
con.execute("create temp table lig(ligand text, inchikey text)"); con.executemany("insert into lig values (?,?)", L.values.tolist())
con.execute("create temp table up(accession text)"); con.executemany("insert into up values (?)", [(u,) for u in ups])
con.execute("create index temp.i1 on lig(inchikey)")
q = """select up.accession, lig.ligand, a.pchembl_value, a.assay_id
from lig join compound_structures cs on cs.standard_inchi_key = lig.inchikey
join activities a on a.molregno = cs.molregno
join assays s on s.assay_id = a.assay_id and s.assay_type = 'B' and s.confidence_score in (8, 9)
join target_dictionary t on t.tid = s.tid and t.target_type = 'SINGLE PROTEIN'
join target_components tc on tc.tid = t.tid
join component_sequences c on c.component_id = tc.component_id
join up on up.accession = c.accession
where a.pchembl_value is not null"""
A = pd.read_sql_query(q, con).rename(columns={"accession": "uniprot"})
P = A.groupby(["uniprot", "ligand"]).agg(pchembl=("pchembl_value", "median"), n_act=("pchembl_value", "size"), n_assays=("assay_id", "nunique")).reset_index()
P.to_csv(f"{DATA}/chembl/ligand_pchembl.csv", index=False)
print(f"{L.inchikey.nunique()} ligand InChIKeys; {len(A)} activities; {len(P)} (UniProt, ligand) pairs with pChEMBL across {P.uniprot.nunique()} UniProts")
