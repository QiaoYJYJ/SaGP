# SaGP
This project builds a classification pipeline using LightGBM and ACC-PSSM features derived from protein sequences.

# Environment Setup
conda create -n pseinone python=2.7

conda activate pseinone

# Step 1: Feature Extraction (ACC-PSSM)
cd Pseinone

python2.7 profile.py ../positiveData.fasta ACC-PSSM -lag 10 -f csv -cpu 64 -out ../feature/ACC-PSSM_pos.csv

python2.7 profile.py ../negativeData.fasta ACC-PSSM -lag 10 -f csv -cpu 64 -out ../feature/ACC-PSSM_neg.csv

This will generate:
feature/ACC-PSSM_pos.csv (positive samples)
feature/ACC-PSSM_neg.csv (negative samples)

# Step 2: Merge Features
Back in the project root directory:

python merge.py (pythhon3 environment)

# Step 3: Prediction
Prepare a feature file for prediction. You must ensure:

The file is named ACC_PSSM.csv

The columns are properly named as ACC_PSSM_F1 to ACC_PSSM_F4000

It contains only the ACC-PSSM features used by the model

Then run:

python predict.py

This script will:

Output results to prediction_results.csv
 
