# PDB `.pkl` → PostgreSQL Pipeline

This repository contains a **step-by-step pipeline** for loading Protein-related data stored in `.pkl` (pickle) files into a **PostgreSQL database**, validating correctness, and exporting/querying the data back into `.pkl` format.

The project is implemented primarily via **Jupyter notebooks**, with detailed explanations and comments embedded directly in the notebooks.

---

## Project Overview

The workflow supports:

- Creating PostgreSQL tables from `.pkl` files
- Exporting complete tables from PostgreSQL back to `.pkl`
- Performing sanity checks to ensure data integrity
- Querying/exporting subsets of data based on `gene_id` and `transcript_id`

This is intended for **data validation, reproducibility, and sharing** of protein-related datasets.

---

## Notebook Descriptions

### 1️⃣ `01_pdb_pkltopsql_pipeline_v02.ipynb`

**Purpose:**  
Core pipeline notebook for database creation and round-trip validation.

**What it does:**
- Loads protein data from a `.pkl` file
- Creates a PostgreSQL table from the `.pkl` data
- Uploads the data into PostgreSQL
- Exports the full table back from PostgreSQL into a `.pkl` file
- Performs basic checks to confirm successful insertion and retrieval

**Notes:**
- This notebook establishes the database schema
- Acts as the foundation for all downstream operations

---

### 2️⃣ `02_pdb_checker_compare_two_pkl_files_v01.ipynb`

**Purpose:**  
Sanity checking and validation of database correctness.

**What it does:**
- Loads the original input `.pkl` file
- Loads the `.pkl` file exported from PostgreSQL
- Compares both datasets
- Verifies that data integrity is preserved across the pkl → PostgreSQL → pkl pipeline

**Notes:**
- Designed to ensure that database operations do not introduce data corruption
- Comparison logic and assumptions are documented in notebook comments

---

### 3️⃣ `03_pdb_filter_geneid_transcriptid_v01.ipynb`

**Purpose:**  
Targeted data extraction from PostgreSQL.

**What it does:**
- Connects to the PostgreSQL database
- Queries data based on:
  - `gene_id`
  - `transcript_id`
- Exports filtered results into `.pkl` files

**Notes:**
- Intended for downstream biological analysis
- Demonstrates selective access patterns on the database

---

## Requirements

The project uses Python with the following core dependencies:

- `pandas`
- `sqlalchemy`
- `psycopg2`
- `jupyterlab` / `notebook`

A full list of dependencies is provided in `requirements.txt`.

---

## Environment Setup (Example)

```bash
conda create -n pdb_env python=3.11
conda activate pdb_env
pip install -r requirements.txt
jupyter lab
