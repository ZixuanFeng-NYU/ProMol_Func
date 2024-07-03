import pandas as pd
import requests
from multiprocessing import Pool

# Function to download CSV for a single AID
def download_csv(aid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?query=download&record_type=datatable&actvty=all&response_type=save&aid={aid}"
    response = requests.get(url)

    if response.status_code == 200:
        with open(f"pubchem_bioassays_with_protein_targets/AID_{aid}.csv", "wb") as f:
            f.write(response.content)
        return f"AID {aid}: CSV file downloaded successfully."
    else:
        return f"AID {aid}: Failed to download CSV file."

if __name__ == "__main__":
    # Read the CSV file
    df = pd.read_csv("pubchem_bioassays_with_protein_target.csv")

    # Get the list of AIDs
    aid_list = df['AID'].tolist()

    # Define the number of processes
    num_processes = 4  # Adjust as needed

    # Create a Pool of processes
    with Pool(num_processes) as pool:
        # Map the download_csv function to each AID in parallel
        results = pool.map(download_csv, aid_list)

    # Print the results
    for result in results:
        print(result)

