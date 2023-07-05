import pandas as pd
import csv
import pymol
# Import PyMOL modules
from pymol import cmd
import os 

#############################################################################################################################

# Making output directory
# Check if the directory already exists
save_dir = os.path.join("Output", "Alignments", "unique_chains")
if not os.path.exists(save_dir):
    os.mkdir(save_dir)

#############################################################################################################################

# Get the directory containing the script
script_directory = os.path.dirname(os.path.abspath(__file__))

file_path = os.path.join(script_directory, "Output", "single_chain_analysis", "unique_chain_data.csv")

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(script_directory) if file.endswith(".pdb")]

# Removing the .pdb extension so it can be matched to the pdb[i]
local_pdb = [filename.replace('.pdb', '') for filename in pdb_files]

pdb = []
chain = []

with open(file_path, 'r') as file:
    # Skip the first line which contains the column name
    next(file)
    reader = csv.reader(file)
    for row in reader:
        pdb.append(row[0])
        chain.append(row[1])


# Set PyMOL to run without GUI
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

# Import First chain of PDBs
for i in range(0, len(pdb)):

    print("Importing " + pdb[i])
    
    if pdb[i] in local_pdb:
        local_pdb_name = pdb[i] + ".pdb"
        cmd.load(local_pdb_name)
    else:
        cmd.fetch(pdb[i], async_=0, quiet=1)

    cmd.remove("{} and not chain {}".format(pdb[i], chain[i]))

    pymol.cmd.save(os.path.join(script_directory, save_dir, "{}.pdb".format(pdb[i] + "_" + chain[i])))

    pymol.cmd.delete(pdb[i])


# Deleting the PDB files that automatically save

# Get the directory containing the script
script_directory = os.path.dirname(os.path.abspath(__file__))

for filename in os.listdir(script_directory):
    if filename.endswith(".cif"):
        file_path = os.path.join(script_directory, filename)
        os.remove(file_path)


# Close PyMol
cmd.quit()
