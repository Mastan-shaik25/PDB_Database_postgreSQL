#!/usr/bin/env python
# coding: utf-8

# # PKL CHECKER TOOLKIT
# - This tool compares input and output .pkl
# - checks if both files are same or not
# - to make sure the code is creating and exporting Proteinn data Tables correctly 

# In[5]:


# =========================================
# PKL CHECKER TOOLKIT
# =========================================
import pickle
import hashlib

# ---------------- CONFIG ---------------- #
INPUT_PKL  = "../data/ghazi/ENSG00000188938.pkl"
OUTPUT_PKL = "../data/ghazi/ENSG00000188938_all_dynamic2.pkl"


# In[6]:


# =========================================
# Utility: Index proteins by (gene_id, transcript_id)
# =========================================
def index_by_key(data):
    """
    Convert a list of proteins to a dict keyed by (gene_id, transcript_id)
    """
    return {(p["gene_id"], p["transcript_id"]): p for p in data}


# In[7]:


# =========================================
# Check 1: Quick sanity (protein counts + keys)
# =========================================
def quick_sanity_check(input_data, output_data):
    print("Number of proteins (input):", len(input_data))
    print("Number of proteins (output):", len(output_data))

    print("Input first keys:", input_data[0].keys())
    print("Output first keys:", output_data[0].keys())
    print("✅ Quick sanity check done\n")

# =========================================
# Check 2: Compare each protein row (nested PDBs)
# =========================================
def compare_proteins(input_data, output_data):
    in_map = index_by_key(input_data)
    out_map = index_by_key(output_data)

    if set(in_map.keys()) != set(out_map.keys()):
        print("❌ Protein key mismatch")
        return False

    for key in in_map:
        inp = in_map[key]
        out = out_map[key]

        # Compare top-level fields
        for field in ["sequence", "exons", "protein_coding", "nmd"]:
            if inp[field] != out[field]:
                print(f"❌ Difference in {field} for {key}")
                return False

        # Compare PDBs
        if len(inp["pdb_files"]) != len(out["pdb_files"]):
            print(f"❌ Different number of PDBs for {key}")
            return False

        for i, (p1, p2) in enumerate(zip(inp["pdb_files"], out["pdb_files"])):
            if p1["pdb_id"] != p2["pdb_id"]:
                print(f"❌ PDB ID mismatch for {key}, index {i}")
                return False
            if p1["content"] != p2["content"]:
                print(f"❌ PDB content mismatch for {key}, index {i}")
                return False

    print("✅ Protein-level check passed")
    return True

# =========================================
# Check 3: Total number of PDBs
# =========================================
def compare_total_pdbs(input_data, output_data):
    input_total = sum(len(p["pdb_files"]) for p in input_data)
    output_total = sum(len(p["pdb_files"]) for p in output_data)
    print("Total PDBs input:", input_total)
    print("Total PDBs output:", output_total)

    if input_total != output_total:
        print("❌ Total PDB count mismatch")
        return False
    print("✅ Total PDB count check passed")
    return True

# =========================================
# Check 4: Full PKL semantic comparison (order-independent)
# =========================================
def compare_pkl_files(file1, file2):
    """
    Compare two PKL files for true semantic identity.
    Checks top-level fields + nested PDBs.
    """
    with open(file1, "rb") as f:
        data1 = pickle.load(f)
    with open(file2, "rb") as f:
        data2 = pickle.load(f)

    map1 = index_by_key(data1)
    map2 = index_by_key(data2)

    if set(map1.keys()) != set(map2.keys()):
        print("❌ Protein key mismatch between files")
        return False

    for key in map1:
        p1 = map1[key]
        p2 = map2[key]

        # Top-level fields
        for field in ["sequence", "exons", "protein_coding", "nmd"]:
            if p1[field] != p2[field]:
                print(f"❌ Difference in '{field}' for protein {key}")
                return False

        # PDBs (order matters)
        if len(p1["pdb_files"]) != len(p2["pdb_files"]):
            print(f"❌ Different number of PDBs for protein {key}")
            return False
        for i, (pb1, pb2) in enumerate(zip(p1["pdb_files"], p2["pdb_files"])):
            if pb1["pdb_id"] != pb2["pdb_id"]:
                print(f"❌ PDB ID mismatch for protein {key}, index {i}")
                return False
            if pb1["content"] != pb2["content"]:
                print(f"❌ PDB content mismatch for protein {key}, index {i}")
                return False

    print("✅ Files match perfectly (semantic + order check)")
    return True

# =========================================
# Check 5: Cryptographic “proof” comparator (byte-level, order-preserving)
# =========================================
def hash_protein(protein):
    h = hashlib.sha256()
    h.update(protein["gene_id"].encode())
    h.update(protein["transcript_id"].encode())
    h.update((protein["sequence"] or "").encode())
    h.update(pickle.dumps(protein["exons"]))
    h.update(str(protein["protein_coding"]).encode())
    h.update(str(protein["nmd"]).encode())
    for pdb in protein["pdb_files"]:
        h.update(pdb["pdb_id"].encode())
        h.update(pdb["content"])
    return h.hexdigest()


def proof_compare_pkl(file1, file2):
    with open(file1, "rb") as f:
        data1 = pickle.load(f)
    with open(file2, "rb") as f:
        data2 = pickle.load(f)

    if len(data1) != len(data2):
        print(f"❌ Different number of proteins: {len(data1)} vs {len(data2)}")
        return False

    for i, (p1, p2) in enumerate(zip(data1, data2)):
        if hash_protein(p1) != hash_protein(p2):
            print(f"❌ Mismatch detected at protein index {i} (order preserved)")
            return False

    print("✅ PROOF PASSED: PKL files are identical (order + content)")
    return True


# In[8]:


# =========================================
# AUTOMATIC EXECUTION (minimal user input)
# =========================================
if __name__ == "__main__":
    # Load data
    with open(INPUT_PKL, "rb") as f:
        input_data = pickle.load(f)
    with open(OUTPUT_PKL, "rb") as f:
        output_data = pickle.load(f)

    # Run all checks
    quick_sanity_check(input_data, output_data)
    compare_total_pdbs(input_data, output_data)
    compare_proteins(input_data, output_data)
    compare_pkl_files(INPUT_PKL, OUTPUT_PKL)
    proof_compare_pkl(INPUT_PKL, OUTPUT_PKL)


# # END OF CHECKS
