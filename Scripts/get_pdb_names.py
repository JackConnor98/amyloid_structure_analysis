import pandas as pd
import csv

df = pd.read_csv('Output/atlas_scraping.txt', delimiter='\t')

# Initialising happy for while loop
try_again = ""

while try_again != "no":
    # Getting user input for pattern matching
    pattern = input("\n Please enter pattern to match your protein of interest: \n\n")

    # filter DataFrame by Protein column matching pattern "syn"
    df_filtered = df[df["Protein"].str.contains(str(pattern), case = False)]

    # print filtered DataFrame
    print(df_filtered)

    # Asking user if they want to try again
    try_again = input("\n Would you like to try again? (yes/no) \n\n")

# Saving filtered DataFrame to a text file
df_filtered.to_csv("Output/selected_pdbs_metadata.txt", index = False, header = True, sep = "\t")

# Extracting PDB names from selected rows
pdb_names = df_filtered["PDB ID"].tolist()

# Printing number of PDBs
print("\n You Selected {} PDBs".format(len(pdb_names)))

# Adding column name to start of the list
pdb_names.insert(0, "PDB")

# Saving PDB names to a text file
with open('Output/pdb_names.txt', 'w') as file:
    # write each value to a new line
    for value in pdb_names:
        file.write(value + '\n')