#  Breast Cancer RNA-Seq Survival Explorer

This interactive web application helps cancer researchers and clinicians explore how gene expression affects survival outcomes in breast cancer subtypes. Built with Dash and Plotly, it enables real-time survival analysis and differential expression comparisons using TCGA and METABRIC datasets — no coding required.

The system has three main components:
1.  A data engineering script (`data_processor.py`) that downloads, cleans, normalizes, and harmonizes RNA-Seq and clinical data from the TCGA and METABRIC cohorts.
2.  An offline deep learning script (`offline_training.py`) that trains a DeepSurv survival model on the TCGA data and uses SHAP (SHapley Additive exPlanations) to interpret the model, generating a ranked list of genes by their prognostic importance.
3.  A real-time backend API (`backend_api.py`) built with Flask and an interactive frontend built with Dash/Plotly. This application serves the pre-computed SHAP rankings and allows users to perform on-the-fly Kaplan-Meier survival analysis and differential expression analysis for any gene, across different molecular subtypes, in both TCGA and METABRIC datasets without requiring any programming knowledge.

##  Features

- **Real-Time Survival Analysis**  
  Select genes and subtypes to generate Kaplan-Meier survival curves with log-rank p-values.

- **Dynamic Differential Expression**  
  Compare tumor subtypes to identify significantly up/down-regulated genes ranked by effect size.

- **Interactive Visualizations**  
  Explore survival curves, box plots, volcano plots, and ranked gene tables with intuitive controls.

- **No Bioinformatics Expertise Needed**  
  Designed for clinicians and researchers with simple dropdowns, sliders, and guided workflows.

## Technologies Used

- [Dash](https://dash.plotly.com/) for web framework and interactivity  
- [Plotly](https://plotly.com/python/) for scientific visualizations  
- [Pandas](https://pandas.pydata.org/) for data manipulation  
- [Gunicorn](https://gunicorn.org/) for deployment

## Data Sources

- **TCGA Breast Cancer Cohort**  
  RNA-Seq and clinical metadata from The Cancer Genome Atlas

- **METABRIC Validation Cohort**  
  Independent dataset for cross-validation and reproducibility

## Installation

```bash
git clone https://github.com/meech253/breast-cancer-dashboard.git
cd breast-cancer-dashboard
pip install -r requirements.txt

python Breast_Cancer_RNA_Seq_Explorer.py
```
## Usage

After starting the server, open http://localhost:8050 in your browser. 

Select a gene and subtype from the dropdowns to generate survival curves and differential expression plots.

## Live Demo
Try the dashboard here: [Breast Cancer RNA‑Seq Explorer](https://breast-cancer-rna-seq-explorer.onrender.com/)

