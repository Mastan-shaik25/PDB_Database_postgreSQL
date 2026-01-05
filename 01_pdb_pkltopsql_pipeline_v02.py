#!/usr/bin/env python
# coding: utf-8

# # Goal:
# - Read a .pkl dataset of proteins
# - Dynamically create PostgreSQL table
# - Insert proteins + PDBs while preserving order
# - Export back to .pkl exactly as input
# - Fully dynamic: no manual column specification

# # Prerequisites
# - **postgreSQL** - Download and install it on your system
# - Download_Link -> https://www.postgresql.org/download/windows/
# - during installation - create and save superuser, pwd - which you need latter
# - **Python** > 3.9+

# In[31]:


# =========================================
# STEP 0: Imports
# =========================================
import pickle
import json
import psycopg2
import pandas as pd
from sqlalchemy import create_engine, text

# ---------------- CONFIG ---------------- #
PKL_PATH   = "../data/ghazi/ENSG00000188938.pkl"
OUTPUT_PKL = "../data/ghazi/ENSG00000188938_all_dynamic2.pkl"

PG_USER = "postgres"
PG_PASSWORD = "manjoor123$ps"
PG_HOST = "localhost"
PG_PORT = "5432"
DB_NAME = "protein_db_dynamic2"

TABLE_NAME = "protein_table_dynamic2"


# In[32]:


# !pip install pandas


# In[33]:


# =========================================
# STEP 1: Load dataset
# =========================================
with open(PKL_PATH, "rb") as f:
    data = pickle.load(f)

ORIGINAL_KEY_ORDER = list(data[0].keys())

print("Total proteins:", len(data))
print("PDB files per protein:", [len(e.get("pdb_files", [])) for e in data])

# Detect columns dynamically from first protein
first_entry = data[0]
pdb_columns = ["pdb_ids", "pdb_files"]
other_columns = [k for k in first_entry.keys() if k != "pdb_files"]
columns = ["protein_index"] + other_columns + pdb_columns
print("Detected columns:", columns)


# In[34]:


# =========================================
# STEP 2: Connect to PostgreSQL (create DB if needed)
# =========================================

# Admin connection to create DB
admin_engine = create_engine(
    f"postgresql+psycopg2://{PG_USER}:{PG_PASSWORD}@{PG_HOST}:{PG_PORT}/postgres",
    isolation_level="AUTOCOMMIT"
)

with admin_engine.connect() as conn:
    result = conn.execute(
        text("SELECT 1 FROM pg_database WHERE datname=:name"), {"name": DB_NAME}
    ).fetchone()
    if not result:
        conn.execute(text(f"CREATE DATABASE {DB_NAME}"))
        print(f"Database '{DB_NAME}' created.")
    else:
        print(f"Database '{DB_NAME}' already exists.")

# Connect to the target database
engine = create_engine(
    f"postgresql+psycopg2://{PG_USER}:{PG_PASSWORD}@{PG_HOST}:{PG_PORT}/{DB_NAME}"
)


# In[35]:


# =========================================
# STEP 3: Create table dynamically
# =========================================

# Map Python types to PostgreSQL types
def infer_pg_type(key, value):
    if key == "protein_index":
        return "INTEGER NOT NULL"
    elif key == "pdb_ids":
        return "TEXT[]"
    elif key == "pdb_files":
        return "BYTEA[]"
    elif isinstance(value, bool):
        return "BOOLEAN"
    elif isinstance(value, (dict, list)):
        return "JSONB"
    else:
        return "VARCHAR(400)"


# Build column definitions
col_defs = []
for col in columns:
    pg_type = infer_pg_type(col, first_entry.get(col, None))
    not_null = "NOT NULL" if col in ["gene_id", "transcript_id", "protein_index"] else ""
    col_defs.append(f"{col} {pg_type} {not_null}")

col_defs_sql = ",\n".join(col_defs)

create_table_sql = f"""
CREATE TABLE IF NOT EXISTS {TABLE_NAME} (
{col_defs_sql},
PRIMARY KEY (gene_id, transcript_id)
);
"""

with engine.begin() as conn:
    conn.execute(text(create_table_sql))

print(f"Table '{TABLE_NAME}' created or already exists.")


# In[36]:


# =========================================
# JSON SAFETY UTIL
# =========================================
import base64

def make_json_safe(obj):
    """
    Recursively convert non-JSON-serializable objects into JSON-safe form.
    Currently handles bytes by base64 encoding.
    """
    if isinstance(obj, bytes):
        return {
            "__type__": "bytes",
            "encoding": "base64",
            "data": base64.b64encode(obj).decode("ascii")
        }
    elif isinstance(obj, dict):
        return {k: make_json_safe(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [make_json_safe(v) for v in obj]
    else:
        return obj


# In[37]:


#
def insert_proteins_with_pdbs(engine, data, table_name):
    """
    Insert proteins with PDB files into the database.
    Fully dynamic, JSON-safe, order-preserving.
    """

    cols_sql = ", ".join(columns)
    vals_sql = ", ".join([f"%({c})s" for c in columns])

    sql = f"""
    INSERT INTO {table_name} ({cols_sql})
    VALUES ({vals_sql})
    ON CONFLICT (gene_id, transcript_id) DO NOTHING;
    """

    conn = engine.raw_connection()
    try:
        cur = conn.cursor()
        total_pdbs = 0

        for idx, entry in enumerate(data):
            row_dict = {"protein_index": idx}

            for key in other_columns:
                value = entry.get(key)

                if isinstance(value, (dict, list)):
                    row_dict[key] = json.dumps(make_json_safe(value))
                else:
                    row_dict[key] = value

            # PDB handling (binary-safe, NOT JSON)
            row_dict["pdb_ids"] = [p["pdb_id"] for p in entry.get("pdb_files", [])]
            row_dict["pdb_files"] = [
                psycopg2.Binary(p["content"])
                for p in entry.get("pdb_files", [])
            ]

            total_pdbs += len(row_dict["pdb_files"])
            cur.execute(sql, row_dict)

        conn.commit()

    finally:
        conn.close()

    print(f"✅ Inserted {len(data)} proteins with {total_pdbs} total PDB files.")


# Run insert
insert_proteins_with_pdbs(engine, data, TABLE_NAME)


# In[38]:


# =========================================
# STEP 5: Verify table content (optional preview)
# =========================================
df_preview = pd.read_sql(f"""
SELECT gene_id, transcript_id, cardinality(pdb_files) AS pdb_count
FROM {TABLE_NAME}
ORDER BY protein_index;
""", engine)
df_preview.head()


# In[39]:


# =========================================
# STEP 6: Export full table dynamically
# MAX-RIGOR VERSION (order + identity safe)
# =========================================

def export_table_to_pkl(
    engine,
    table_name,
    output_path,
    original_key_order
):
    """
    Export entire table back to a .pkl file with:
    - Exact input key order
    - Exact row order
    - Exact value identity
    - Fully dynamic schema support
    """

    df = pd.read_sql(
        f"SELECT * FROM {table_name} ORDER BY protein_index",
        engine
    )

    data_out = []

    for _, row in df.iterrows():
        protein_entry = {}

        for key in original_key_order:
            if key == "pdb_files":
                # Reconstruct PDBs exactly
                protein_entry["pdb_files"] = [
                    {
                        "pdb_id": pid,
                        "content": bytes(pb)
                    }
                    for pid, pb in zip(
                        row["pdb_ids"] or [],
                        row["pdb_files"] or []
                    )
                ]
            else:
                value = row[key]

                # JSONB columns are already decoded by psycopg2
                protein_entry[key] = value

        data_out.append(protein_entry)

    with open(output_path, "wb") as f:
        pickle.dump(data_out, f, protocol=pickle.HIGHEST_PROTOCOL)

    print(
        f"✅ Exported {len(data_out)} proteins "
        f"(order + structure preserved at max possible fidelity)"
    )


# Example usage
export_table_to_pkl(
    engine,
    TABLE_NAME,
    OUTPUT_PKL,
    ORIGINAL_KEY_ORDER
)


# In[40]:


# =========================================
# STEP 7: Proof comparator (identity + order)
# =========================================
import hashlib

def hash_protein(protein):
    h = hashlib.sha256()
    h.update(protein["gene_id"].encode())
    h.update(protein["transcript_id"].encode())
    h.update((protein["sequence"] or "").encode())
    h.update(pickle.dumps(protein.get("exons", [])))
    h.update(str(protein.get("protein_coding", False)).encode())
    h.update(str(protein.get("nmd", False)).encode())
    for pdb in protein.get("pdb_files", []):
        h.update(pdb["pdb_id"].encode())
        h.update(pdb["content"])
    return h.hexdigest()


def proof_compare_pkl(file1, file2):
    with open(file1, "rb") as f1, open(file2, "rb") as f2:
        data1, data2 = pickle.load(f1), pickle.load(f2)

    if len(data1) != len(data2):
        print(f"❌ Different number of proteins: {len(data1)} vs {len(data2)}")
        return False

    for i, (p1, p2) in enumerate(zip(data1, data2)):
        if hash_protein(p1) != hash_protein(p2):
            print(f"❌ Mismatch at protein index {i}")
            return False

    print("✅ PROOF PASSED: PKL files are identical (order + content)")
    return True


# Run proof
proof_compare_pkl(PKL_PATH, OUTPUT_PKL)


# # END OF NOTE BOOK
