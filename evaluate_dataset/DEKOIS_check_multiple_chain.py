import os

# Directory containing the FASTA files
input_dir = "DEKOIS2.0_library_protein"

# Get the list of FASTA files
list_of_fasta = os.listdir(input_dir)

# Iterate through each FASTA file
for file_ in list_of_fasta:
    if file_.endswith(".fasta"):
        fasta_path = os.path.join(input_dir, file_)

        with open(fasta_path, "r") as f:
            # Read all lines in the file
            lines = f.readlines()

        # Variables to track the current sequence and the longest sequence
        current_header = None
        current_sequence = []
        longest_header = None
        longest_sequence = ""

        for line in lines:
            if line.startswith('>'):
                # Save the previous sequence if it was longer than the current longest
                if current_header and len(''.join(current_sequence)) > len(longest_sequence):
                    longest_header = current_header
                    longest_sequence = ''.join(current_sequence)
                
                # Start a new sequence
                current_header = line.strip()
                current_sequence = []
            else:
                # Add line to the current sequence
                current_sequence.append(line.strip())

        # Final check for the last sequence in the file
        if current_header and len(''.join(current_sequence)) > len(longest_sequence):
            longest_header = current_header
            longest_sequence = ''.join(current_sequence)

        # Write the longest sequence back to the file
        with open(fasta_path, "w") as f_out:
            if longest_header and longest_sequence:
                f_out.write(f"{longest_header}\n")
                f_out.write(f"{longest_sequence}\n")

        # Print the result
        print(f"Processed {file_}: longest sequence length {len(longest_sequence)}")

