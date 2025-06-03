import pandas as pd
from Bio import SeqIO

# === Step 1: Load FASTA IDs ===
def extract_ids(fasta_path):
    return [record.id for record in SeqIO.parse(fasta_path, "fasta")]

neg_ids = extract_ids("../negativeData.fasta")
pos_ids = extract_ids("../positiveData.fasta")

# === Step 2: Load feature CSVs ===
neg_df = pd.read_csv("ACC-PSSM_neg.csv", header=None)
pos_df = pd.read_csv("ACC-PSSM_pos.csv", header=None)

# === Step 3: Add ID and Label columns ===
neg_df.insert(0, "ID", neg_ids)
neg_df["Label"] = 0

pos_df.insert(0, "ID", pos_ids)
pos_df["Label"] = 1

# === Step 4: Concatenate ===
merged_df = pd.concat([neg_df, pos_df], ignore_index=True)

# === Step 5: Generate header names ===
num_features = merged_df.shape[1] - 2  # exclude ID and Label
feature_cols = [f"ACC_PSSM_F{i+1}" for i in range(num_features)]
merged_df.columns = ["ID"] + feature_cols + ["Label"]

# === Step 6: Save to CSV ===
merged_df.to_csv("filtered_features.csv", index=False)
print("✅ Merged and saved as filtered_features.csv")