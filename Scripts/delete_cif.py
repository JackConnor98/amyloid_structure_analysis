import os

# Deleting the PDB files that automatically save

# Get the directory containing the script
script_directory = os.path.dirname(os.path.abspath(__file__))

# Set the directory to the script directory
directory = script_directory

for filename in os.listdir(directory):
    if filename.endswith(".cif"):
        file_path = os.path.join(directory, filename)
        os.remove(file_path)