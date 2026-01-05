#!/usr/bin/env python
# coding: utf-8

# # Protein Filter Export Notebook
# # This notebook allows you to:
# - Filter proteins by `gene_id` or `gene_id + transcript_id`
# - Export complete rows (including all PDB files) as a `.pkl` file
# - Fully standalone: just update `gene_id` and `transcript_id`  

# In[15]:


import pickle
import json
import pandas as pd
import psycopg2
from sqlalchemy import create_engine, text


# In[16]:


# ================= USER CONFIG =================
#
GENE_ID = "ENSG00000188938"         # mention geneid to export
TRANSCRIPT_ID = "ENST00000375412"   # optional
OUTPUT_DIR = "../data/ghazi/"

PG_USER = "postgres"
PG_PASSWORD = "manjoor123$ps"
PG_HOST = "localhost"
PG_PORT = "5432"
DB_NAME = "protein_db_dynamic1"
TABLE_NAME = "protein_table_dynamic1"
# ==============================================


# In[17]:


# Connect to Database
engine = create_engine(
    f"postgresql+psycopg2://{PG_USER}:{PG_PASSWORD}@{PG_HOST}:{PG_PORT}/{DB_NAME}"
)


# In[18]:


# Detect Original PKL Schema Dynamically
#
def detect_other_columns(engine, table_name):
    """
    Detect original PKL columns dynamically.
    Excludes DB-only fields.
    """
    df = pd.read_sql(
        f"SELECT * FROM {table_name} ORDER BY protein_index LIMIT 1",
        engine
    )

    exclude = {"protein_index", "pdb_ids", "pdb_files"}
    return [c for c in df.columns if c not in exclude]


OTHER_COLUMNS = detect_other_columns(engine, TABLE_NAME)
print("Detected PKL columns:", OTHER_COLUMNS)


# In[19]:


# Canonical Export Helper - This is the same logic as the pipeline exporter
#
def rows_to_pkl(df, other_columns, output_path):
    data_out = []

    for _, row in df.iterrows():
        # Safety check
        assert len(row["pdb_ids"] or []) == len(row["pdb_files"] or [])

        protein_entry = {}

        for col in other_columns:
            if col == "exons":
                protein_entry[col] = row[col] or []
            else:
                protein_entry[col] = row[col]

        protein_entry["pdb_files"] = [
            {"pdb_id": pid, "content": bytes(pb)}
            for pid, pb in zip(row["pdb_ids"] or [], row["pdb_files"] or [])
        ]

        data_out.append(protein_entry)

    with open(output_path, "wb") as f:
        pickle.dump(data_out, f)

    print(f"✅ Exported {len(data_out)} proteins → {output_path}")


# In[20]:


# Filter Function 1 — by gene_id
#
def export_by_gene(engine, table_name, gene_id, other_columns, output_path):
    df = pd.read_sql(
        f"""
        SELECT * FROM {table_name}
        WHERE gene_id = %s
        ORDER BY protein_index
        """,
        engine,
        params=(gene_id,)
    )

    if df.empty:
        print(f"❌ No rows found for gene_id={gene_id}")
        return

    rows_to_pkl(df, other_columns, output_path)

# Filter Function 2 — by gene_id + transcript_id
#
def export_by_gene_transcript(engine, table_name, gene_id, transcript_id, other_columns, output_path):
    df = pd.read_sql(
        f"""
        SELECT * FROM {table_name}
        WHERE gene_id = %s AND transcript_id = %s
        ORDER BY protein_index
        """,
        engine,
        params=(gene_id, transcript_id)
    )

    if df.empty:
        print(f"❌ No row found for gene_id={gene_id}, transcript_id={transcript_id}")
        return

    rows_to_pkl(df, other_columns, output_path)


# In[21]:


# Run Exports
# Export all proteins for a gene
export_by_gene(
    engine,
    TABLE_NAME,
    GENE_ID,
    OTHER_COLUMNS,
    f"{OUTPUT_DIR}/{GENE_ID}_gene.pkl"
)

# Export single transcript (optional)
export_by_gene_transcript(
    engine,
    TABLE_NAME,
    GENE_ID,
    TRANSCRIPT_ID,
    OTHER_COLUMNS,
    f"{OUTPUT_DIR}/{GENE_ID}_{TRANSCRIPT_ID}.pkl"
)


# # END OF CHECKS
