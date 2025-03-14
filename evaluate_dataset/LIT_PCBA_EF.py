import pandas as pd
import os
import math


# Get the list of files in the directory
listfile = os.listdir("LIT_PCBA_EF/data_per_target")
Targets,results = [],[]  # Store results for saving
for file in listfile:
    if file.endswith("ProMol_Func_prediction.csv"):
        target = file.split("/")[-1].split("_")[0]
        filepath = os.path.join("LIT_PCBA_EF/data_per_target", file)

        # Load the CSV file
        df = pd.read_csv(filepath)

        # Sort by predictions in descending order (higher is better)
        df = df.sort_values(by="preds", ascending=False).reset_index(drop=True)

        # Total number of molecules
        total_mols = len(df)
        total_positives = df['Targets'].sum()

        # Define percentages
        percentages = [0.005,0.01, 0.05, 0.10]

        # Calculate enrichment factors
        enrichment_factors = {}
        for pct in percentages:
            top_n = math.ceil(total_mols * pct)   # Ensure at least one sample
            top_hits = df.iloc[:top_n]['Targets'].sum()
            enrichment_factors[f"EF_{pct*100:.1f}%"] = (top_hits / total_positives) / pct if total_positives > 0 else 0

        print(f"Target: {target}, Enrichment Factors: {enrichment_factors}")
        # Append results
        Targets.append(target)
        results.append(enrichment_factors)

# Save results to CSV
if results:
    df_results = pd.DataFrame(results)
    df_results.insert(0, "Target", Targets)  # Add 'Target' as a column instead of index
    df_results.to_csv("LIT_PCBA_EF/enrichment_factors.csv", index=False)
    print("Enrichment factors saved to 'LIT_PCBA_EF/LIT_PCBA_EF_per_target.csv'")
~                                                                             
