import pandas as pd
import numpy as np
import os
import sys
from scipy import stats
from statsmodels.stats.multitest import multipletests

# --- PATH CONFIGURATION ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")

def find_column_insensitive(df, target_names):
    """Helper to find a column looking for multiple name variations."""
    for col in df.columns:
        if col in target_names:
            return col
        # Case-insensitive check
        if col.upper() in [t.upper() for t in target_names]:
            return col
    return None

def load_and_clean_tcga():
    print("--- Loading TCGA Data ---")
    
    # 1. Load Clinical Data
    clinical_path = os.path.join(DATA_DIR, "brca_tcga_pan_can_atlas_2018", "data_clinical_patient.txt")
    if not os.path.exists(clinical_path):
        print(f"ERROR: Could not find file at {clinical_path}")
        return None

    clinical = pd.read_csv(clinical_path, sep="\t", comment="#", header=0)
    
    # Clean Survival Data
    clinical['status'] = clinical['OS_STATUS'].apply(lambda x: 1 if 'DECEASED' in str(x) else 0)
    clinical['time'] = pd.to_numeric(clinical['OS_MONTHS'], errors='coerce')
    
    # 2. Handle Subtypes (The robust part)
    # Try to find subtype in patient file first
    subtype_col = find_column_insensitive(clinical, ['Subtype', 'SUBTYPE', 'Subtype_Selected'])
    
    if subtype_col:
        print(f"Found subtype column in patient file: {subtype_col}")
        clinical = clinical.rename(columns={subtype_col: 'subtype'})
    else:
        # Fallback: Load sample file
        print("Subtype not in patient file. Checking sample file...")
        sample_path = os.path.join(DATA_DIR, "brca_tcga_pan_can_atlas_2018", "data_clinical_sample.txt")
        
        if os.path.exists(sample_path):
            sample_df = pd.read_csv(sample_path, sep="\t", comment="#", header=0)
            
            # Look for subtype in sample file
            sample_subtype_col = find_column_insensitive(sample_df, ['Subtype', 'SUBTYPE', 'Subtype_Selected'])
            
            if sample_subtype_col:
                print(f"Found subtype column in sample file: {sample_subtype_col}")
                # Merge it into clinical
                sample_df = sample_df.rename(columns={sample_subtype_col: 'subtype'})
                clinical = clinical.merge(sample_df[['PATIENT_ID', 'subtype']], on='PATIENT_ID', how='left')
            else:
                print("WARNING: Could not find 'Subtype' column even in sample file.")
                print(f"Available columns: {sample_df.columns.tolist()}")
                return None # Stop here if we can't find subtypes
        else:
             print(f"WARNING: Sample file not found at {sample_path}")

    # Check if 'subtype' exists now
    if 'subtype' not in clinical.columns:
        print("CRITICAL ERROR: 'subtype' column missing after all attempts.")
        return None

    clinical = clinical[['PATIENT_ID', 'time', 'status', 'subtype']].dropna(subset=['time', 'status'])
    clinical.set_index('PATIENT_ID', inplace=True)

    # 3. Load Expression Data
    # Try the raw file first, then fall back to z-scores version
    expr_path = os.path.join(DATA_DIR, "brca_tcga_pan_can_atlas_2018", "data_mrna_seq_v2_rsem.txt")
    if not os.path.exists(expr_path):
        # Use z-scores version if raw file doesn't exist
        expr_path = os.path.join(DATA_DIR, "brca_tcga_pan_can_atlas_2018", "data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt")
        if not os.path.exists(expr_path):
            print(f"ERROR: Could not find expression file at {expr_path}")
            return None
        else:
            print(f"Using z-scores expression file: {expr_path}")
         
    expression = pd.read_csv(expr_path, sep="\t", index_col='Hugo_Symbol')
    expression = expression.loc[~expression.index.duplicated(keep='first')]
    expression = expression.drop(columns=['Entrez_Gene_Id'], errors='ignore').T
    
    # Fix IDs
    expression.index = expression.index.str[:-3]
    
    # 4. Merge
    master_df = clinical.join(expression, how='inner')
    print(f"TCGA SUCCESS: Loaded {master_df.shape[0]} patients, {master_df.shape[1]} columns")
    return master_df

# --- EXECUTION BLOCK ---
if __name__ == "__main__":
    tcga_df = load_and_clean_tcga()
    
    if tcga_df is not None:
        print("\nPreview of TCGA Data:")
        print(tcga_df[['time', 'status', 'subtype']].head())
        
        if 'TP53' in tcga_df.columns:
            print(f"\nTP53 Gene Found! Mean Expression: {tcga_df['TP53'].mean():.4f}")