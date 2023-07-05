import pandas as pd
import csv
import pymol
# Import PyMOL modules
from pymol import cmd
import os 

#############################################################################################################################

# Making output directory
# Check if the directory already exists
save_dir = os.path.join("Output", "Clustering", "grouped_alignments")
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

#############################################################################################################################

# Get the directory containing the script
script_directory = os.path.dirname(os.path.abspath(__file__))

# Get the path to the .csv file containing the cluster information
file_path = os.path.join(script_directory, "Output", "Clustering", "grouped_pdb_ids.csv")

# Get the path to the folder containing all the unique PDBs
unique_chain_folder = os.path.join(script_directory, "Output", "Alignments", "unique_chains")

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(script_directory) if file.endswith(".pdb")]

# Removing the .pdb extension so it can be matched to the pdb[i]
local_pdb = [filename.replace('.pdb', '') for filename in pdb_files]

# Initializing lists
pdb = []
chain = []
group = []

with open(file_path, 'r') as file:
    # Skip the first line which contains the column name
    next(file)
    reader = csv.reader(file)
    for row in reader:
        pdb.append(row[0])
        chain.append(row[1])
        group.append(row[3])


# Set PyMOL to run without GUI
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()


# Setting import to one group at a time
for i in range(1, int(max(group))+1):
    
    print("Group", i)

    # Making group directory
    group_dir = os.path.join(save_dir, "group_{}".format(i))
    if not os.path.exists(group_dir):
        os.mkdir(group_dir)

    # Looping through all PDBs
    for j in range(0, len(pdb)):

        # Only import chains that match the current group
        if group[j] == str(i): 

            pdb_id = pdb[j] + "_" + chain[j]

            print("Importing " + pdb_id)
            
            chain_path = os.path.join(unique_chain_folder, pdb_id + ".pdb")

            if pdb[i] in local_pdb:
             local_pdb_name = pdb[i] + ".pdb"
         
            cmd.load(chain_path, pdb_id)

    # Aligning chains
    pymol.cmd.extra_fit(selection="(all)", reference=None, method="align")

    # Saving Alignment as a PDB
    pymol.cmd.save(os.path.join(group_dir, "group_{}_alignment.pse".format(i))) 

    # Clearing environment for next group
    pymol.cmd.delete("all")

# Close PyMol
cmd.quit()
