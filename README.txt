================================================================
**README.txt for Breast Cancer RNA-Seq Explorer (Team 205)**
================================================================

DESCRIPTION
-----------
This package contains the complete source code and supporting data for the "Breast Cancer RNA-Seq Explorer," a project for CSE6242 / CX4242. The project consists of an AI-driven pipeline and an interactive web application designed to identify prognostic biomarkers in breast cancer and allow for real-time survival analysis.

The system has three main components:
1.  A data engineering script (`data_processor.py`) that downloads, cleans, normalizes, and harmonizes RNA-Seq and clinical data from the TCGA and METABRIC cohorts.
2.  An offline deep learning script (`offline_training.py`) that trains a DeepSurv survival model on the TCGA data and uses SHAP (SHapley Additive exPlanations) to interpret the model, generating a ranked list of genes by their prognostic importance.
3.  A real-time backend API (`backend_api.py`) built with Flask and an interactive frontend built with Dash/Plotly. This application serves the pre-computed SHAP rankings and allows users to perform on-the-fly Kaplan-Meier survival analysis and differential expression analysis for any gene, across different molecular subtypes, in both TCGA and METABRIC datasets without requiring any programming knowledge.


INSTALLATION
------------
The project is built in Python 3.9+. To set up the environment and install all necessary dependencies, please follow these steps:

1.  **Clone the repository or unzip the CODE folder.**
    Navigate to the `CODE/` directory in your terminal.

2.  **Create a virtual environment (recommended):**
    python -m venv venv
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`

3.  **Install dependencies from requirements.txt:**
    A `requirements.txt` file is provided with all necessary packages. Run the following command:
    pip install -r requirements.txt

    This will install all required libraries, including: pandas, numpy, torch, torchtuples, pycox, shap, lifelines, scikit-learn, Flask, Dash, and Plotly.

4.  **Download Data (handled by the script):**
    The data processing script is configured to download the necessary data from cBioPortal automatically. No manual data download is required. The script will create a `data/` directory to store these files.


EXECUTION
---------
There are three main scripts that can be executed.

**1. To Preprocess Data and Train the AI Model:**
   This step runs the complete data processing and offline training pipeline. It will download data, clean it, train the DeepSurv model, and generate the `gene_importance.csv` file.

   Open a terminal in the `CODE/` directory and run:
   python offline_training.py

   Expected output: The script will print progress messages for data loading, preprocessing, model training, and SHAP value calculation. Upon completion, a file named `gene_importance.csv` will be created in the `CODE/` directory. This step can take several minutes depending on your hardware.

**2. To Run the Interactive Web Application (Main Demo):**
   This is the primary execution step to launch the interactive application. This script relies on the output from the training step (`gene_importance.csv`).

   Open a terminal in the `CODE/` directory and run:
   python backend_api.py

   Expected output: The terminal will show that the Flask/Dash server is running and will provide a local URL, typically http://127.0.0.1:8050/.

   To use the application, open a web browser (like Chrome or Firefox) and navigate to:
   [**http://127.0.0.1:8050/**](http://127.0.0.1:8050/)

   You will see the interactive dashboard where you can explore gene expression, perform survival analysis, and view the AI-driven gene rankings.

**3. To Run Only the Data Processing Step:**
   If you only wish to run the data preprocessing part without training the model, you can modify the `offline_training.py` script by commenting out the function calls related to model training. However, the recommended execution path is to run `offline_training.py` followed by `backend_api.py`.

# Github link
https://github.com/meech253/breast-cancer-dashboard

# Dashboard link
https://breast-cancer-rna-seq-explorer.onrender.com
