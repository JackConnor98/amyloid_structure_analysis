
# Pymol is installed elsewhere so brings up error but code runs fine
import pymol
# Import PyMOL modules
from pymol import cmd
import os 

#############################################################################################################################

# Making output directory
# Check if the directory already exists
if not os.path.exists("Output"):
    os.mkdir("Output")

if not os.path.exists("Output/Alignments"):
    os.mkdir("Output/Alignments")

#############################################################################################################################

# Reading in all the .cif files in the directory

# Get the path to the current directory
folder_path = os.path.dirname(os.path.realpath(__file__))

# Get all the.cif files
cif_files = [file for file in os.listdir(folder_path) if file.endswith(".cif")]

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(folder_path) if file.endswith(".pdb")]

all_files = cif_files + pdb_files

# Set PyMOL to run without GUI
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()


for i in all_files:

    # Setting name of loaded structures to be the pdb code
    object_name = i.split(".")[0]

    # Loading in the structure and setting the name 
    cmd.load(i, object_name)

    # Get the list of chains
    chain_list = cmd.get_chains(object_name)

    # Unload structures 
    cmd.remove(object_name)

    # Seperating structures into individual chains
    for chain in chain_list:

        # Setting name of loaded structures to be the pdb code + the chain 
        chain_name = object_name + "_" + chain
        print(chain_name)

        # Re-loading the structure and setting the name 
        cmd.load(i, chain_name)

        # Removing all but one chain
        cmd.remove("{} and not chain {}".format(chain_name, chain))

        # Getting the list of residues in the chain
        atoms = cmd.get_model(chain_name).atom
        
        # Check if chain contains unkown residues
        unknown = False
        for atom in atoms:
            if atom.resn == "UNK":
                unknown = True
                break

        # Unload chains that contain unkown residues
        if unknown == True:
            cmd.remove(chain_name)

# Aligning all single chains from the same pdb file
# Aligning everything
pymol.cmd.extra_fit(selection='(all)', reference=None, method='align')

# Saving Alignment
pymol.cmd.save("Output/Alignments/all_PDBs_all_chains_aligned.cif")


# Close PyMol
cmd.quit()
