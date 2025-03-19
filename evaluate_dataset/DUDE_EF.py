import pandas as pd
import os
import math

# Get the list of files in the directory
listfile = os.listdir("DUDE_per_target_data")
Targets, results = [], []  # Store results for saving

for file in listfile:
    if file.endswith("ProMol_Func_prediction.csv"):
        target = file.split("/")[-1].split("_")[0]
        filepath = os.path.join("DUDE_per_target_data", file)

        # Load the CSV file
        df = pd.read_csv(filepath)

        # Sort by predictions in descending order (higher is better)
        df = df.sort_values(by="preds", ascending=False).reset_index(drop=True)

        # Total number of molecules
        total_mols = len(df)
        total_positives = df['Targets'].sum()

        # Define percentages
        percentages = [0.005, 0.01, 0.05, 0.10]

        # Calculate enrichment factors
        enrichment_factors = {}
        for pct in percentages:
            top_n = math.ceil(total_mols * pct)  # Ensure at least one sample
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

    # Calculate mean enrichment factors across all targets
    mean_values = df_results.iloc[:, 1:].mean().to_dict()
    mean_values["Target"] = "Mean"  # Label the mean row

    # Append mean row to the dataframe
    df_results = pd.concat([df_results, pd.DataFrame([mean_values])], ignore_index=True)

    df_results.to_csv("DUDE_per_target_data/DUDE_EF_results.csv", index=False)
    print("Enrichment factors saved to 'DUDE_per_target_data/DUDE_EF_results.csv' with mean values included.")
