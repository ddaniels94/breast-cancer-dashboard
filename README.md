#  Breast Cancer RNA-Seq Survival Explorer

This interactive web application helps cancer researchers and clinicians explore how gene expression affects survival outcomes in breast cancer subtypes. Built with Dash and Plotly, it enables real-time survival analysis and differential expression comparisons using TCGA and METABRIC datasets — no coding required.

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

python Breast
