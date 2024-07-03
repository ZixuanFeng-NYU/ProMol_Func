import requests
import pandas as pd
from multiprocessing import Pool

def fetch_data(i):
    try:
        print(i)
        aid = i
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/assay/" + str(i) + "/JSON"
        response = requests.get(url)

        if response.status_code == 200:
            data = response.json()

            current_assay_type = "N/A"
            current_target = "N/A"
            current_target_url = "N/A"

            for item in data.get('Record', {}).get('Section', []):
                if item.get('TOCHeading') == 'Target':
                    current_target = item['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
                    try:
                        current_target_url = item['Section'][0]['Information'][0]['URL']
                    except:
                        current_target_url = 'N/A'
                elif item.get('TOCHeading') == 'BioAssay Annotations':
                    current_assay_type = item['Information'][0]['Value']['StringWithMarkup'][0]['String']
            #print("AID:", aid, "assay_type:", current_assay_type, "Target:", current_target, "Target_URL:", current_target_url)
            return {'AID': aid, 'assay_type': current_assay_type, 'Target': current_target, 'Target_URL': current_target_url}
        else:
            print("Error:", response.status_code)
            return {'AID': aid, 'assay_type': "N/A", 'Target': "N/A", 'Target_URL': "N/A"}

    except Exception as e:
        print(f"Error processing AID {i}: {str(e)}")
        return {'AID': aid, 'assay_type': "N/A", 'Target': "N/A", 'Target_URL': "N/A"}

if __name__ == "__main__":
    with Pool() as pool:
        results = pool.map(fetch_data, range(1, 1964001))

    # Create a DataFrame to store the collected data
    df = pd.DataFrame(results)

    # Save the DataFrame to a CSV file
    df.to_csv("pubchem_assay_type.csv", index=False)

