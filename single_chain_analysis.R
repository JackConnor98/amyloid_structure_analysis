cat("Running single_chain_analysis.R\n")

# Installing and Loading Packages

#install.packages(c("dplyr", "stringr", "tidyr", "ggplot2", "scales", "readr"))

library("dplyr")
library("stringr")
library("tidyr")
library("ggplot2")
library("scales")
library("readr")

############################################################################################

# Creating Directories
if(!dir.exists("Output/single_chain_analysis")) dir.create("Output/single_chain_analysis")

############################################################################################

# Read in df
df <- read.csv("Output/chain_alignment_data.csv")

# Ordering by residue_position 
df <- df %>% arrange(residue_position)


# Removing non-fibril core resiudes

metadata <- read_delim("Output/selected_pdbs_metadata.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

metadata <- subset(metadata, select = c("PDB ID",  "Residues Ordered"))
names(metadata) <- c("pdb", "ordered_residues")

metadata$core_start <- NA
metadata$core_end <- NA

for (i in 1:nrow(metadata)) {
  
  range <- metadata$ordered_residues[i]
  range <- str_split_1(range, paste("-", ",", sep = "|"))
  
  metadata$core_start[i] <- as.numeric(range[1])
  metadata$core_end[i] <- as.numeric(range[length(range)])
    
}

# Combining df and metadata
df <- merge(df, metadata, all = TRUE)

# Changing Start and End NAs for non-amyloid atlas pdbs
df <- df %>%
  mutate(core_start = replace_na(core_start, 0),
         core_end = replace_na(core_end, Inf))

# Filtering residues that are outside the core
df <- df %>% filter(residue_position >= core_start, 
                    residue_position <= core_end)


# Assigning Molecular Weights to Each Atom
molecular_weights <- data.frame(atom_id = c("C", "H", "O", "N", "S"), mw = c(12, 1, 16, 14, 32))

# Adding MW to df
df <- merge(df, molecular_weights, by = "atom_id")

# Cleaning Environment
rm(list = setdiff(ls(), c("df")))
invisible(gc())

##############################################################################################

# Getting amyloid names
amyloid_names <- unique(df$pdb)

###############################################################################################
###############################################################################################
###############################################################################################

# End to end distance analysis

min_pos <- df %>% 
            group_by(pdb, chain) %>% 
              slice_min(residue_position, n = 1) %>% # slice_min() returns entire row of min residue_position
                select(pdb, chain, residue_position, x, y, z) # Keeps columns of interest

names(min_pos)[names(min_pos) == "residue_position"] <- "first_residue"
names(min_pos)[names(min_pos) == "x"] <- "first_x"
names(min_pos)[names(min_pos) == "y"] <- "first_y"
names(min_pos)[names(min_pos) == "z"] <- "first_z"

max_pos <- df %>% 
            group_by(pdb, chain) %>% 
             slice_max(residue_position, n = 1) %>% # slice_min() returns entire row of min residue_position
              select(pdb, chain, residue_position, x, y, z) # Keeps columns of interest

names(max_pos)[names(max_pos) == "residue_position"] <- "last_residue"
names(max_pos)[names(max_pos) == "x"] <- "last_x"
names(max_pos)[names(max_pos) == "y"] <- "last_y"
names(max_pos)[names(max_pos) == "z"] <- "last_z"

e2e_data <- merge(min_pos, max_pos, by = c("pdb", "chain"))

e2e_data <- e2e_data %>%
              mutate(e2e_distance = sqrt((first_x - last_x)^2 + (first_y - last_y)^2 + (first_z - last_z)^2))

# Cleaning Environment
rm(list = c("min_pos", "max_pos"))
invisible(gc())

###############################################################################################
###############################################################################################
###############################################################################################

# Radius of Gyration Formula
# RoG = RMSD from the centre of mass 

# Centre of Mass Calculation (CoM)

#x_CoM = 1/total_chain_mass * sum(mass*x_coord)

# Multiplying each atom's mw by its coordinate
df <- df %>% group_by(pdb, chain) %>% mutate(mass_x_x = mw * x,
                                             mass_x_y = mw * y,
                                             mass_x_z = mw * z)

# calculating the sum of mass * coordinate for each chain
sum_x_coord <- df %>% group_by(pdb, chain) %>% summarise_at(vars(mass_x_x), funs(sum))
sum_y_coord <- df %>% group_by(pdb, chain) %>% summarise_at(vars(mass_x_y), funs(sum))
sum_z_coord <- df %>% group_by(pdb, chain) %>% summarise_at(vars(mass_x_z), funs(sum))

# Finding the total mass of each chain
total_chain_mass <- df %>% group_by(pdb, chain) %>% summarise_at(vars(mw), funs(sum))
names(total_chain_mass)[names(total_chain_mass) == "mw"] <- "total_chain_mass"

# Creating a Centre of Mass data frame
CoM_data <- merge(merge(merge(sum_x_coord, sum_y_coord, by = c("pdb", "chain"), all = TRUE),
                              sum_z_coord, by = c("pdb", "chain"), all = TRUE),
                              total_chain_mass, by = c("pdb", "chain"), all = TRUE)

CoM_data <- CoM_data %>% mutate(x_CoM = (1/total_chain_mass) * mass_x_x,
                                y_CoM = (1/total_chain_mass) * mass_x_y,
                                z_CoM = (1/total_chain_mass) * mass_x_z)

# Cleaning Environment
rm(list = c("sum_x_coord", "sum_y_coord", "sum_z_coord", "total_chain_mass"))

# Creating Radius of Gyration
Rg_data <- merge(df, CoM_data, by = c("pdb", "chain"), all = TRUE)

# Calculating the distance of each atom from the CoM
Rg_data <- Rg_data %>% mutate(distance_to_CoM = sqrt((x_CoM - x)^2 + (y_CoM - y)^2 + (z_CoM - z)^2))

Rg_data <- Rg_data %>% mutate(mass_x_squared_distance_to_CoM = (mw * (distance_to_CoM)^2))

# Calculating the sum of mass * distance^2 for each atom
Rg_data <- Rg_data %>% group_by(pdb, chain) %>% summarise_at(vars(mass_x_squared_distance_to_CoM), funs(sum))

# Adding total chain mass and CoM Coordinates
Rg_data$x_CoM <- CoM_data$x_CoM
Rg_data$y_CoM <- CoM_data$y_CoM
Rg_data$z_CoM <- CoM_data$z_CoM
Rg_data$total_chain_mass <- CoM_data$total_chain_mass

# Calculating Rg by dividing the sum of mass * distance^2 by the total mw
Rg_data <- Rg_data %>% group_by(pdb, chain) %>% mutate(Rg = sqrt(mass_x_squared_distance_to_CoM / total_chain_mass))

# Dropping mass_x_squared_distance_to_CoM
Rg_data <- subset(Rg_data, select = -c(mass_x_squared_distance_to_CoM))

###############################################################################################

# Creating output_data used for plotting and saving

output_data <- Rg_data
output_data$e2e_distance <- e2e_data$e2e_distance

# Cleaning Environment
rm(list = c("CoM_data", "e2e_data", "Rg_data", "mean_data"))

###############################################################################################

# Plotting 

# Getting the last 3 column names for plotting
vars_to_plot <- c(names(output_data)[(ncol(output_data)-2):ncol(output_data)])

# Creating names for axis labels
axis_label <- c("Total Chain Mass", "Radius of Gyration", "End to End Distance")


for (i in 1:length(vars_to_plot)) {
  
  ggplot(output_data, aes(x = pdb, y = .data[[vars_to_plot[i]]], colour = pdb)) + # .data - refers to the current data frame being used
    geom_point(size = 4) +                                                        # [[]] is used to extract the column
    ylab(paste(axis_label[i])) +
    xlab("PDB") +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
          panel.grid.major = element_line(colour = "grey", linewidth = 0.8, linetype = "dashed"),
          panel.grid.minor = element_line(colour = "grey", linewidth = 0.8, linetype = "dashed"),
          plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
          axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
          axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle = 45, hjust=1),
          axis.title=element_text(size = 25, face = "bold"),
          axis.line.x = element_line(color="black", linewidth = 1.5),
          axis.line.y = element_line(color="black", linewidth = 1.5),
          legend.position = "none",
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 10, face = "bold")) 
  
  ggsave(plot = last_plot(), file = paste0("Output/single_chain_analysis/", vars_to_plot[i], ".png"), 
        height = 10, width = 20)
  
  
}

###############################################################################################
###############################################################################################
###############################################################################################

# Finding Outliers

# Finding mean for each pdb
mean_data <- output_data %>% 
                group_by(pdb) %>% 
                  summarise_at(vars(total_chain_mass, Rg, e2e_distance), funs(mean))

# Rename Columns
names(mean_data) <- c("pdb", "mean_total_chain_mass", "mean_Rg", "mean_e2e_distance")

# Adding means back to output_data 
output_data <- merge(output_data, mean_data, by = "pdb") 

output_data <- output_data %>% 
                 mutate(norm_total_chain_mass = total_chain_mass / mean_total_chain_mass,
                        norm_Rg = Rg / mean_Rg,
                        norm_e2e_distance = e2e_distance / mean_e2e_distance)

# Saving output_data
write.csv(output_data, "Output/single_chain_analysis/single_chain_analysis.csv", row.names = FALSE)

##############################################################################################

# Normalised Plotting

# Getting the last 3 column names for plotting
vars_to_plot <- c(names(output_data)[(ncol(output_data)-2):ncol(output_data)])

# Creating names for axis labels
axis_label <- c("Total Chain Mass", "Radius of Gyration", "End to End Distance")

# Define the color gradient
colours <- c("blue", "black", "red")  # Example colors for the gradient

for (i in 1:length(vars_to_plot)) {
  
  ggplot(output_data, aes(x = pdb, y = .data[[vars_to_plot[i]]], colour = .data[[vars_to_plot[i]]])) + # .data - refers to the current data frame being used
    geom_hline(yintercept = c(0.5,0.9,1,1.1,1.5), 
               colour = c("red", "orange", "black", "orange", "red"),
               linetype = c("dashed", "dashed", "solid", "dashed", "dashed"),
               linewidth = 3) +
    geom_point(size = 6) +                                  
    ylab(paste("Normalised ", axis_label[i])) +
    xlab("PDB") +
    scale_y_continuous(breaks = seq(0,10,0.25)) +
    scale_color_gradientn(colours = colours, values = rescale(c(0, 0.5, 1))) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
          panel.grid.major = element_line(colour = "grey", linewidth = 0.8, linetype = "dashed"),
          panel.grid.minor = element_line(colour = "grey", linewidth = 0.8, linetype = "dashed"),
          plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
          axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
          axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle = 45, hjust=1),
          axis.title=element_text(size = 25, face = "bold"),
          axis.line.x = element_line(color="black", linewidth = 1.5),
          axis.line.y = element_line(color="black", linewidth = 1.5),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 12, face = "bold")) +
    labs(colour = paste(axis_label[i]))
  
  ggsave(plot = last_plot(), file = paste0("Output/single_chain_analysis/", vars_to_plot[i], ".png"), 
         height = 10, width = 20)
  
  
}

# Cleaning Environment
rm(list = c("axis_label", "colours", "vars_to_plot", "i"))
invisible(gc())

##############################################################################################
##############################################################################################
##############################################################################################

#######################################
# Identifying all the distinct chains #
#######################################

# Select values that have over 20% difference between min and max

# Selecting norm data columns from output_data
norm_data <- output_data[, c(1:2, (ncol(output_data) - 2):ncol(output_data))]

# Finding the minimum and maximum norm values for each pdb
diff_max_min <- norm_data %>% group_by(pdb) %>% summarise_at(c("norm_e2e_distance", "norm_Rg", "norm_total_chain_mass"),
                                                    funs(min, max))

# Finding the difference between Max and Min norm values
diff_max_min <- diff_max_min %>% mutate(diff_e2e = round(norm_e2e_distance_max - norm_e2e_distance_min, 2),
                                        diff_Rg = round(norm_Rg_max - norm_Rg_min, 2),
                                        diff_total_chain_mass = round(norm_total_chain_mass_max - norm_total_chain_mass_min, 2))

# Selecting columns of interest
diff_max_min <- diff_max_min[, c(1, (ncol(diff_max_min) - 2):ncol(diff_max_min))]

# Selecting only PDBs where max-min difference is greater than 20% for any metric 
#variable_pdbs <- diff_max_min %>% filter(if_sum(c(diff_e2e, diff_Rg, diff_total_chain_mass), ~ . > 0.2))

# Selecting only PDBs where max-min difference is greter than 20% for 2+ metrics 
variable_pdbs <- diff_max_min %>%
                  filter((diff_e2e > 0.2) + (diff_Rg > 0.2) + (diff_total_chain_mass > 0.2) > 1)


##############################################################################################

#######################################
# Creating minimum unique chains data #
#######################################

# Getting variable pdb names
variable_pdb_names <- variable_pdbs[,1]

# Merging them to norm_data to get back all the chain names and norm values
merged <- merge(norm_data, variable_pdb_names, by = "pdb")

# Rounding norm_data to 2dp
merged$norm_e2e_distance <- round(merged$norm_e2e_distance, 1)
merged$norm_Rg <- round(merged$norm_Rg, 1)
merged$norm_total_chain_mass <- round(merged$norm_total_chain_mass, 1)

# Removing all but one identical chains, grouped by PDB (not including e2e_dist as its too sensitive)
variable_chains <- merged %>% 
                    group_by(pdb) %>%
                     distinct(norm_total_chain_mass, norm_Rg, .keep_all = TRUE)

variable_chain_names <- variable_chains[,1:2]

# Cleaning Environment
rm(list = c("diff_max_min", "mean_data", "merged", "norm_data", "variable_pdbs", "variable_pdb_names", "variable_chains"))
invisible(gc())

# Now we have a list of each PDBs with multiple distinct chains and one of each distinct chain
# Next we need to filter df to remove all except chain A for non variable PDBs + the identified distinct chains

# Removing previous calculation columns
df <- df[1:(ncol(df)-4)]

# Creating a df containing only the variable chains
variable_df <- merge(variable_chain_names, df, by = c("pdb", "chain"))

# Identifying chain names from each PDB
chain_names <- df %>% group_by(pdb) %>% distinct(chain)
chain_names <- chain_names %>% arrange(pdb, chain)
single_chain <- distinct(chain_names, pdb, .keep_all = TRUE)

# Selecting only the first chain for each PDB from df 
constant_df <- df 

for (i in 1:nrow(single_chain)) {
  
  constant_df <- constant_df %>% filter(!(pdb==single_chain$pdb[i] & chain != single_chain$chain[i]))
  
}

# Move pdb and chain to the start of min_df
#constant_df <- constant_df %>% select(pdb, chain, everything()) # Doesn't work in VSCode

# Removing the variable PDBs, leaving only the PDBs with no significant variation
constant_df <- constant_df %>% filter(!(pdb %in% unique(variable_df$pdb)))

# Creating min_df with all the atom coordinates for one chain from each PDB + any identified variable chains
#min_df <- rbind(constant_df, variable_df) # rbind() doesn't run in VSCode
min_df <- merge(constant_df, variable_df, by = colnames(constant_df), all = TRUE)

# Order min_df by PDB name, chain and residue_position
min_df <- min_df[order(min_df$pdb, min_df$chain, min_df$residue_position),]

# Creating pdb_id so if PDB has multiple distinct chains they are saved as PDB_chain
min_df$pdb_id <- ""

for (i in 1:nrow(min_df)) {
  
  if (min_df$chain[i] == single_chain$chain[single_chain$pdb == min_df$pdb[i]]) {
    
    min_df$pdb_id[i] <- as.character(min_df$pdb[i]) # As character is needed for it to work in VSCode, not needed in RStudio

  } else {
    
    min_df$pdb_id[i] <- paste0(as.character(min_df$pdb[i]), "_", as.character(min_df$chain[i]))
    
  }
  
}

# Saving min_df
write.csv(min_df, file = "Output/unique_chains_alignment.csv", row.names = FALSE)

# Saving only the PDB name and chain
unique_pdbs <- subset(min_df, select = c("pdb", "chain"))
unique_pdbs <- subset(unique_pdbs, !duplicated(min_df[c("pdb", "chain")]))
write.csv(unique_pdbs, file = "Output/unique_chains.csv", row.names = FALSE)


# Saving Unique chains E2E, Rg and total_mws

# Extracting unique pdb and chain names
unique_pdb_chains <- subset(min_df, select = c("pdb", "chain"))

# Extracting data for unique chains
unique_data <- merge(unique_pdb_chains, output_data, by = c("pdb", "chain"))
unique_data <- distinct(unique_data)

#  Selecting columns of interest
unique_data <- subset(unique_data, select = c("pdb", "chain", "e2e_distance", "Rg", "total_chain_mass"))

# Saving dataframe
write.csv(unique_data, file = "Output/single_chain_analysis/unique_chain_data.csv", row.names = FALSE)


# Cleaning Environment
rm(list = c("variable_chain_names", "variable_df", "constant_df", "unique_pdb_chains", "unique_pdbs", "i"))
invisible(gc())

###############################################################################################

# Getting .pdb names for RMSD Calculations
pdb_names <- list.files(pattern = "\\.pdb$")
pdb_names <- gsub("\\.pdb$", "", pdb_names)

###############################################################################################

# Plotting E2E Distance

ggplot(unique_data, aes(x = reorder(pdb, e2e_distance), y = e2e_distance)) +
  geom_point(aes(color = pdb %in% pdb_names), size = 4) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
  ggtitle("End-to-End Distance") +
  ylab("E2E Distance") +
  xlab("PDB") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(unique_data$e2e_distance) * 1.05)) +
  theme(panel.grid.major = element_line(colour = "grey30", linewidth = 0.5, 
                                        linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5), 
        axis.title = element_text(size = 18, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black",
                                   angle = 90, vjust = 0.4, hjust = 1),
        axis.text.y = element_text(size = 14, colour = "black"))

ggsave(plot = last_plot(), file = "Output/single_chain_analysis/e2e_plot.png", width = 14, height = 8)


###############################################################################################

# Plotting total chain mass

ggplot(unique_data, aes(x = reorder(pdb, total_chain_mass), y = total_chain_mass)) +
  geom_point(aes(color = pdb %in% pdb_names), size = 4) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
  ggtitle("Total Chain Mass") +
  ylab("Total Chain Mass") +
  xlab("PDB") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(unique_data$total_chain_mass) * 1.05)) +
  theme(panel.grid.major = element_line(colour = "grey30", linewidth = 0.5, 
                                        linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5), 
        axis.title = element_text(size = 18, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black",
                                   angle = 90, vjust = 0.4, hjust = 1),
        axis.text.y = element_text(size = 14, colour = "black"))

ggsave(plot = last_plot(), file = "Output/single_chain_analysis/total_mw_plot.png", width = 14, height = 8)


###############################################################################################


# Plotting Rg data

ggplot(unique_data, aes(x = reorder(pdb, Rg), y = Rg)) +
  geom_point(aes(color = pdb %in% pdb_names), size = 4) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
  ggtitle("Radius of Gyration") +
  ylab("Rg") +
  xlab("PDB") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(unique_data$Rg) * 1.05)) +
  theme(panel.grid.major = element_line(colour = "grey30", linewidth = 0.5, 
                                        linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5), 
        axis.title = element_text(size = 18, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black",
                                   angle = 90, vjust = 0.4, hjust = 1),
        axis.text.y = element_text(size = 14, colour = "black"))

ggsave(plot = last_plot(), file = "Output/single_chain_analysis/Rg_plot.png", width = 14, height = 8)


