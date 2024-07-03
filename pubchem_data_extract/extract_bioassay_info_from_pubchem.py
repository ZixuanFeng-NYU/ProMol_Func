import requests
import pandas as pd

aid, assay_type, target, target_url = [], [], [], []

for i in range(1, 1964001):
    print(i)
    aid.append(i)
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/assay/" + str(i) + "/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()

        # Initialize placeholders for current assay's information
        current_assay_type = "N/A"
        current_target = "N/A"
        current_target_url = "N/A"

        # Iterate over sections of the JSON response
        for item in data.get('Record', {}).get('Section', []):
            if item.get('TOCHeading') == 'Target':
                print(item)
                current_target = item['Section'][0]['Information'][0]['Value']['StringWithMarkup'][0]['String']
                print(current_target)
                try:
                    current_target_url = item['Section'][0]['Information'][0]['URL']
                    print(current_target_url)
                except:
                    current_target_url='N/A'
                    print(current_target_url)
            elif item.get('TOCHeading') == 'BioAssay Annotations':
                print(item)
                current_assay_type = item['Information'][0]['Value']['StringWithMarkup'][0]['String']
                print(current_assay_type)
        # Append the collected information to respective lists
        assay_type.append(current_assay_type)
        target.append(current_target)
        target_url.append(current_target_url)

    else:
        print("Error:", response.status_code)
        # Append placeholders for failed requests
        assay_type.append("N/A")
        target.append("N/A")
        target_url.append("N/A")

# Create a DataFrame to store the collected data
df = pd.DataFrame({
    'AID': aid,
    'assay_type': assay_type,
    'Target': target,
    'Target_URL': target_url
})

# Save the DataFrame to a CSV file
df.to_csv("pubchem_assay_type.csv", index=False)

