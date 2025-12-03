import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split as tts
from sklearn.metrics import mean_squared_error
from sklearn.ensemble import RandomForestRegressor

class Project_X:
    def load_and_align_data(self, gene_file_path, drug_file_path):
        print("Loading Data... this might take a minute.")
        try:
            # 1. Load Gene Expression
            # Using index_col=0 to handle the gene names
            gene_df = pd.read_csv(gene_file_path, sep='\t', index_col=0)
            
            # Clean columns: "DATA.906826" -> "906826"
            gene_df.columns = gene_df.columns.str.replace('DATA.', '')
            
            # 2. Load Drug Response
            drug_raw = pd.read_csv(drug_file_path)
            
            # Pivot to Matrix: Rows=Cells, Cols=Drugs
            drug_df = drug_raw.pivot_table(index='Cosmic ID', columns='Drug Name', values='IC50', aggfunc='mean')
            
            # 3. Alignment
            # Ensure indices are strings for matching
            gene_df.columns = gene_df.columns.astype(str)
            drug_df.index = drug_df.index.astype(str)
            
            # Find intersection
            common_cells = drug_df.index.intersection(gene_df.columns)
            print(f"Found {len(common_cells)} common cell lines.")
            
            # Filter and Transpose Gene Data (So Rows=Cells, Cols=Genes)
            X = gene_df[common_cells].T
            
            # Filter Drug Data (So Rows=Cells, Cols=Drugs)
            Y = drug_df.loc[common_cells]
            
            print(f"Final X Shape (Genetics): {X.shape}")
            print(f"Final Y Shape (Drugs):    {Y.shape}")
            
            return X, Y

        except FileNotFoundError as e:
            print(e)
            return None, None
        
            
    def find_biomarkers(self, X, Y, drug_name):
        try:
            if drug_name not in Y.columns:
                    print(f"Error: Drug '{drug_name}' not found in dataset.")
                    return None
            
            drug_series = Y[drug_name]

            # Handle missing data(NaN)
            drug_series.dropna()
            valid_cells = drug_series.dropna().index

            y_clean = drug_series.loc[valid_cells]
            x_clean = X.loc[valid_cells]

            # Calculate correlation
            correlation_matrix = x_clean.corrwith(y_clean)

            # Select top 50 genes by absolute value
            abs_corr_matrix = correlation_matrix.abs()
            sorted_corr_matrix = abs_corr_matrix.sort_values(ascending = False)
            top_50_genes = sorted_corr_matrix.head(50).index.tolist()

            print(f"Top 5 genes found: {top_50_genes[: 5]}")

            return top_50_genes, correlation_matrix
        
        except Exception as e:
            print(f"An error occurred during correlation analysis: {e}")
            return None, None
        
    def train_model(self, X, Y, top_genes, drug_name):
        try:
            # Handle Drug Name missing error
            if drug_name not in Y.columns:
               print(f"Error: {drug_name} not in drug matrix.")
               return None
            
            # Taking all of the genes will be exhaustive, so we consider only genes whose expresstion is most affected
            X_subset = X[top_genes]
            Y_target = Y[drug_name]

            # Handling missing values(NaN)
            valid_indices = Y_target.dropna().index

            # Selecting only the valid values and removing the rest
            X_clean = X_subset.loc[valid_indices]
            Y_clean = Y_target.loc[valid_indices]

            X_train, X_test, Y_train, Y_test = tts(X_clean, Y_clean, test_size=0.2, random_state=42)

            Model = RandomForestRegressor(n_estimators=100, random_state=42)

            Model.fit(X_train, Y_train)

            # Calculate the r2 value (r2 value > 0.3 is acceptable)
            R_squared_value = Model.score(X_test, Y_test)

            # Calculate the root mean squared error:
            Y_predicted = Model.predict(X_test)
            rmse = np.sqrt(mean_squared_error(Y_test, Y_predicted))

            # Print the results
            print(f"Model Performance for {drug_name}:")
            print(f"  R-squared: {R_squared_value:.4f}")
            print(f"  RMSE: {rmse:.4f}")

            return Model, R_squared_value, rmse

        except Exception as e:
           print(f"An unexpected error occured:{e}")
           return None, None, None

sol = Project_X()

# Following files should be downloaded from the official GDSC website (https://www.cancerrxgene.org/downloads/anova)
# as the files are too large to upload here
gene_file = 'Cell_line_RMA_proc_basalExp.txt'
drug_file = 'PANCANCER_IC_Tue Nov 25 07_52_02 2025.csv'

X, Y = sol.load_and_align_data(gene_file_path=gene_file, drug_file_path=drug_file)

if X is not None:
    target_drug = "Cisplatin"
    top_genes, _ = sol.find_biomarkers(X, Y, target_drug)
            
    if top_genes:
        # Train the model for 'Cisplatin'
        sol.train_model(X, Y, top_genes, target_drug)
else:
    print("Error: One or both data files are missing.") # Missing file error handling


# ------------- Thank you GDSC for the data files ----------------