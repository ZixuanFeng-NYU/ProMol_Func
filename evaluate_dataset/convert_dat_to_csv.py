import csv

# Function to convert .dat file to .csv file
def dat_to_csv(dat_file_path, csv_file_path, delimiter=' '):
    with open(dat_file_path, 'r') as dat_file:
        dat_reader = csv.reader(dat_file, delimiter=delimiter)
        
        with open(csv_file_path, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file)
            for row in dat_reader:
                csv_writer.writerow(row)

# Example usage
dat_file_path = 'CASF-2016/power_screening/CoreSet.dat'  # Replace with your .dat file path
csv_file_path = 'CASF-2016/power_screening/CoreSet.csv'  # Replace with your desired .csv file path
dat_to_csv(dat_file_path, csv_file_path, delimiter=' ')  # Replace delimiter if needed (e.g., ',' for comma, '\t' for tab)

