cat("Running create_df.R\n")

# Installing and Loading Packages
library("dplyr")
library("stringr")
library("tidyr")
library("ggplot2")
library("data.table")

############################################################################################

# Creating Directories
#if(!dir.exists("Output/Data")) dir.create("Output/Data")
#if(!dir.exists("Output/Figures")) dir.create("Output/Figures")

############################################################################################

# Importing Data

cat("\nReading in Data\n")

# Read all the lines of a file
lines <- readLines("Output/Alignments/all_PDBs_all_chains_aligned.cif")

# Initialize empty vectors for amyloid names and atom lines
amyloid_names <- c()
atom_lines <- c()

# Initialize new variables to hold the current name and row before adding them to the list
current_name <- ""
current_row <- ""

# Progress update
progress <- 0

# Loops through lines of .cif file
# If line begins with data, extract the name (name of pdb file)
# every line beginning with ATOM until the next data line will contain the coordinates for the pdb file
# add the data name to these rows so we have a column specifying which pdb the coordinates belong to
# upon hitting a new data line we want to change the current_name to the next pdb file name etc
for (line in lines) {
  if (substr(line, 1, 4) == "data") {
    
    current_name <- ""
    
    # Removing data_ from name
    current_name <- gsub("data_", "", line)
    
    # Removing _X from end of name (where X = chain)
    current_name <- substr(current_name, 1, nchar(current_name) - 2)
    
    # Adding current name to the names list
    amyloid_names <- append(amyloid_names, current_name)
  }
  
  if (substr(line, 1, 4) == "ATOM") {
    
    # Removing ATOM from the start of each row
    current_row <- gsub("ATOM   ", "", line)
    
    # Separating row into columns by space
    current_row <- paste(current_row, current_name, sep=" ")
    
    # Adding current row to list
    atom_lines <- append(atom_lines, current_row)
  }

  # Displaying progress of reading in data
  #progress <- progress + 1
  #percent <- progress / length(lines) * 100
  #percent <- round(percent, 2)

  #if (is.integer(percent) == TRUE) {
  #  cat("Progress: ", percent, " %\n")
  #}
  
  #cat("Progress: ", percent, " %\n")
  
}

# Cleaning Environment
rm(list = c("current_name", "current_row", "line", "lines"))
invisible(gc())


# Since each row in atom_lines uses different amount of spaces as separators 
# Split lines into fields based on one or more spaces
atom_lines_list <- strsplit(atom_lines, "\\s+")

# Convert list of fields into a data frame
cat("Creating Dataframe\n")
df <- as.data.frame(do.call(rbind, atom_lines_list), stringsAsFactors = FALSE)
# Creates a matrix where column 1 will contain the first element from the list with each row being a different list (i.e. line)
# do.call(rbind, atom_lines_list)

############################################################################################

# data.table's fread() is very fast but inst as flexible and can't create pdb column
#df <- fread("Output/Alignments/all_PDBs_all_chains_aligned.cif")

# Data Tidying
cat("Tidying Data")

# Find column names that contain only "." or "?"
cols_to_delete <- colnames(df)[apply(df == "." | df == "?", 2, all)] # apply(df == "." | df == "?", 2, all): This applies the all function column-wise (2) to the logical matrix generated in the previous step, resulting in a logical vector that indicates whether all values in each column are equal to "." or "?".

# Remove the columns from the data frame
df <- df[, -which(colnames(df) %in% cols_to_delete)] # does this by finding the indices of the column whose 
#names match those in cols_to_delete

# Removing unwanted columns
df<- df[, -which(colnames(df) %in% c("V4", "V6", "V7"))]


# Renaming Columns
df_col_names <- c("atom_number", "atom_id", "atom_alt_conformation", "residue_name", "residue_position", 
                  "x", "y", "z", "occupancy", "temp_factor", "charge", "chain", "segment_identifier", "pdb")
names(df) <- df_col_names

df$atom_number <- as.numeric(df$atom_number)
df$residue_position <- as.numeric(df$residue_position)
df$x <- as.numeric(df$x)
df$y <- as.numeric(df$y)
df$z <- as.numeric(df$z)

# Cleaning Environment
rm(list = c("atom_lines", "atom_lines_list", "cols_to_delete", "df_col_names"))


# Filtering out all positions except central Carbon
df <- df %>% filter(atom_alt_conformation == "CA")

# Removing rows where a pdb file has two alpha-C's for a single residue position
df <- df %>% distinct(residue_position, chain, pdb, .keep_all = TRUE)

# Saving df
write.csv(df, "Output/chain_alignment_data.csv", row.names = FALSE)
