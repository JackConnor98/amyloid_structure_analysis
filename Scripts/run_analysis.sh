
#!/bin/bash

# Setting Run Parameters
scrape=1                       # 0 - Don't Web Scrape | 1 - Web Scrape Amyloid Atlas
pdb_selection=1                 # 0 - Don't Select | 1 - Select
align=1                         # 0 - Don't Align | 1 - Do Alignment
delete_cif=1                    # 0 - Keep .CIF files | 1 - Delete .CIF files
analysis=1                      # 0 - Do not run analysis | 1 - Run analysis

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

if [ $scrape -eq  1 ]; then
    # Web scraping Amyloid Atlas
    python amyloid_atlas_scraper.py
fi

if [ $pdb_selection -eq 1 ]; then
    # Select which PDB files to run
    python get_pdb_names.py
fi

if [ $align -eq 1 ]; then
    # Fetch pdbs to get .cif files
    python fetch_pdb.py

    # Aligning all chains from all pdb files
    python PDB_alignment.py

    if [ $delete_cif -eq 1 ]; then
        # Remove all CIF files
        python delete_cif.py
    fi

    # Creating the alignment dataframe
    Rscript create_df.R
fi

if [ $analysis -eq 1 ]; then

    # Creating Ordered Residues Plots
    Rscript plot_ordered_residues.R

    # Calculate the Rg and E2E_distance for all chains
    Rscript single_chain_analysis.R

    # Saving unique chains
    python save_unique_chains.py

    # Calculate the RMSD between all unique PDB structures
    Rscript RMSD_analysis.R

    # Clustering Analysis
    Rscript clustering.R

    # Aligning unique chains by cluster group
    python group_cluster_alignment.py
    
fi
