import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats
from statsmodels.stats.multitest import multipletests
from data_processor import load_and_clean_tcga
import os

# --- INITIALIZATION ---
print("Initializing Backend Engine...")
TCGA_DF = load_and_clean_tcga()

# Load AI Importance Scores
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
IMPORTANCE_FILE = os.path.join(BASE_DIR, "gene_importance.csv")
try:
    IMPORTANCE_DF = pd.read_csv(IMPORTANCE_FILE).set_index("Gene")
    print("AI Importance Scores Loaded.")
except FileNotFoundError:
    print(f"WARNING: gene_importance.csv not found at {IMPORTANCE_FILE}. Run offline_training.py first.")
    IMPORTANCE_DF = pd.DataFrame()

# Export IMPORTANCE_DF for use in dashboard
__all__ = ['get_survival_analysis', 'get_gene_importance', 'get_expression_by_subtype', 'TCGA_DF', 'SUBTYPE_MAP', 'IMPORTANCE_DF']

# --- SUBTYPE TRANSLATION LAYER (The Fix) ---
# This maps human-readable names to the internal data codes
SUBTYPE_MAP = {
    'TCGA': {
        'Basal': 'BRCA_Basal',
        'Basal-like': 'BRCA_Basal',
        'Luminal A': 'BRCA_LumA',
        'LumA': 'BRCA_LumA',
        'Luminal B': 'BRCA_LumB',
        'LumB': 'BRCA_LumB',
        'HER2': 'BRCA_Her2',
        'HER2-enriched': 'BRCA_Her2',
        'Normal': 'BRCA_Normal'
    },
}

def get_survival_analysis(gene_name, target_subtype, dataset='TCGA'):
    """
    Generates Kaplan-Meier data for a specific gene and subtype.
    Uses TCGA dataset.
    """
    # Select Dataset
    df = TCGA_DF
    
    if df is None: return {"error": "Dataset not loaded"}
    
    # Translate Subtype Name (e.g. "Basal-like" -> "BRCA_Basal")
    dataset_map = SUBTYPE_MAP.get(dataset, {})
    clean_subtype = dataset_map.get(target_subtype, target_subtype) # Default to input if not found
    
    # Filter by Subtype
    subset = df[df['subtype'] == clean_subtype].copy()
    
    if subset.empty:
        # Debugging help: Show what subtypes ARE available
        available = df['subtype'].unique().tolist()
        return {"error": f"No samples for '{target_subtype}' (Mapped to: '{clean_subtype}'). Available: {available}"}
        
    if gene_name not in subset.columns: 
        return {"error": f"Gene {gene_name} not found in {dataset}"}
    
    # Stratify: High vs Low
    median = subset[gene_name].median()
    high_mask = subset[gene_name] >= median
    
    T_high = subset.loc[high_mask, 'time']
    E_high = subset.loc[high_mask, 'status']
    T_low = subset.loc[~high_mask, 'time']
    E_low = subset.loc[~high_mask, 'status']
    
    # Log-Rank Test
    try:
        p_value = logrank_test(T_high, T_low, event_observed_A=E_high, event_observed_B=E_low).p_value
    except:
        p_value = 1.0 # If fails (e.g., no events), assume no diff
        
    # KM Curves
    kmf = KaplanMeierFitter()
    kmf.fit(T_high, E_high)
    high_curve = kmf.survival_function_
    
    kmf.fit(T_low, E_low)
    low_curve = kmf.survival_function_
    
    return {
        "gene": gene_name,
        "subtype": target_subtype,
        "mapped_subtype": clean_subtype,
        "p_value": p_value,
        "plot_data": {
            "high": {"time": high_curve.index.tolist(), "prob": high_curve.iloc[:,0].tolist()},
            "low": {"time": low_curve.index.tolist(), "prob": low_curve.iloc[:,0].tolist()}
        }
    }

def get_gene_importance(gene_name):
    """
    Returns the SHAP score from your Offline AI model.
    """
    if not IMPORTANCE_DF.empty and gene_name in IMPORTANCE_DF.index:
        score = IMPORTANCE_DF.loc[gene_name, "SHAP_Importance"]
        return {"gene": gene_name, "shap_score": round(float(score), 5)}
    return {"gene": gene_name, "shap_score": "N/A"}

def get_expression_by_subtype(gene_name, dataset='TCGA'):
    """
    Returns expression data for a gene broken down by subtype for box plot visualization.
    """
    df = TCGA_DF
    
    if df is None:
        return {"error": "Dataset not loaded"}
    
    if gene_name not in df.columns:
        return {"error": f"Gene {gene_name} not found in {dataset}"}
    
    # Get all subtypes and their expression values
    expression_data = {}
    reverse_map = {}
    
    # Build reverse map for human-readable names
    if dataset in SUBTYPE_MAP:
        for human_name, internal_code in SUBTYPE_MAP[dataset].items():
            reverse_map[internal_code] = human_name
    
    # Group by subtype
    for subtype_code in df['subtype'].dropna().unique():
        subset = df[df['subtype'] == subtype_code]
        expression_values = subset[gene_name].dropna().tolist()
        
        if len(expression_values) > 0:
            # Use human-readable name if available
            display_name = reverse_map.get(subtype_code, str(subtype_code))
            expression_data[display_name] = expression_values
    
    return {
        "gene": gene_name,
        "expression_data": expression_data
    }

# --- Test Block ---
if __name__ == "__main__":
    print("\n--- Testing API ---")
    # Test 1: Survival (We use 'Basal' now, which maps to BRCA_Basal automatically)
    print(get_survival_analysis("TP53", "Basal", dataset='TCGA'))
    
    # Test 2: AI Score
    print(get_gene_importance("COL1A1"))