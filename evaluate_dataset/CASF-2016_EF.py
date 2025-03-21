import pandas as pd
import numpy as np
import math

# Load the dataset
df = pd.read_csv("CASF-2016_organized_data_ProMol_Func_prediction.csv")  # Update with the actual file path

# Function to calculate EF for each pro_id
def calculate_ef(df, percentages=[0.005, 0.01, 0.05, 0.10]):
    ef_results = []

    # Process each pro_id separately
    for pro_id, group in df.groupby('pro_id'):
        group = group.sort_values(by='preds', ascending=False)  # Sort by prediction scores

        total_compounds = len(group)  # Total compounds for this pro_id
        total_actives = (group['Targets'] == 1).sum()  # Total actives
        total_decoys = (group['Targets'] == 0).sum()  # Total decoys

        ef_dict = {
            'pro_id': pro_id,
            'total_actives': total_actives,
            'total_decoys': total_decoys,
            'total_compounds': total_compounds
        }

        # Calculate EF for each percentage
        for pct in percentages:
            num_selected = max(1, math.ceil(total_compounds * pct))  # Ensure at least 1 compound
            top_selected = group.head(num_selected)  # Get top x% compounds
            actives_in_top = (top_selected['Targets'] == 1).sum()  # Count actives in top x%

            # Compute EF (handle case when total_actives is zero)
            ef = ((actives_in_top / num_selected) / (total_actives / total_compounds)) if total_actives > 0 else 0
            ef_dict[f'EF{int(pct*100)}%'] = round(ef, 4)  # Round for better readability

        ef_results.append(ef_dict)

    ef_df = pd.DataFrame(ef_results)

    # Add a final row with mean EF values
    mean_values = ef_df.mean(numeric_only=True).to_dict()
    mean_values['pro_id'] = 'Mean'  # Label the mean row
    ef_df = pd.concat([ef_df, pd.DataFrame([mean_values])], ignore_index=True)

    return ef_df

# Compute EF for each pro_id
ef_df = calculate_ef(df)

# Save results
ef_df.to_csv("CASF-2016_EF_results_with_counts.csv", index=False)

# Print results
print(ef_df)
