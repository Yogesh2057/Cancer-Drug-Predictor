# Cancer-Drug-Predictor

Cancer Drug Response Predictor (Machine Learning & Systems Biology)

This project implements a full-stack Machine Learning pipeline to predict cancer cell line sensitivity to specific drugs based on their baseline gene expression profiles.

By integrating genomics (Gene Expression) and pharmacogenomics (Drug Response) data, the model identifies molecular biomarkers and trains a regressor to predict IC50 values (drug concentration required to kill 50% of cells).

🧬 Biological Validation: The "Camptothecin" Case Study

To validate the pipeline, the model was trained on Camptothecin (a Topoisomerase I inhibitor).

Model Performance: Achieved an R-squared of 0.49 and RMSE of 1.37, indicating strong predictive power for a biological system.

Biomarker Discovery: The model independently identified SLFN11 (Schlafen 11) as the top predictor of sensitivity.

Scientific Context: High expression of SLFN11 is a clinically validated biomarker for sensitivity to DNA-damaging agents (like Camptothecin), confirming that the model is capturing real biological signals, not just statistical noise.

🛠 Technical Architecture

The project is structured into three modular components:

1. Data Harvesting & Engineering (ETL)

Source: Genomics of Drug Sensitivity in Cancer (GDSC) Database.

Processing:

Parses massive gene expression matrices (.txt, ~17,000 genes) and drug response tables (.csv).

Implements Pivot Tables to transform long-format drug data into matrices.

Aligns genomic and pharmacological datasets using index intersection to ensure perfect sample matching (~900 common cell lines).

2. Feature Selection (Biomarker Discovery)

Dimensionality Reduction: Instead of using all 17,000 genes, the script performs a Pearson Correlation Analysis between the drug's IC50 vector and every gene expression vector.

Selection: Automatically selects the top 50 genes with the highest absolute correlation (Sensitivity or Resistance biomarkers) to prevent overfitting.

3. Machine Learning (Predictive Modeling)

Algorithm: Random Forest Regressor (scikit-learn).

Workflow:

Aligns X (Top 50 Genes) and Y (Drug IC50), handling missing values dynamically.

Splits data into Training (80%) and Testing (20%) sets.

Evaluates performance using RMSE and R-squared metrics on unseen test data.

📂 Data Sources

Gene Expression: RMA-Normalized Basal Expression (GDSC).

Drug Response: GDSC2 Fitted Dose Response (IC50).

Note: While the code logic is original, the raw biological data is sourced from the public GDSC repository (Wellcome Sanger Institute).

🤖 Acknowledgements

Project Architecture & Guidance: Developed with the assistance of Google Gemini, which provided the roadmap for the system architecture, debugging support, and theoretical explanations of the machine learning workflows.
