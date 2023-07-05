cat("Running RMSD_analysis.R\n")

# Installing and Loading Packages
library("dplyr")
library("stringr")
library("tidyr")
library("ggplot2")

############################################################################################

# Creating Directories
if(!dir.exists("Output/interPDB_analysis")) dir.create("Output/interPDB_analysis")
if(!dir.exists("Output/interPDB_analysis/RMSD")) dir.create("Output/interPDB_analysis/RMSD")
if(!dir.exists("Output/interPDB_analysis/Distance")) dir.create("Output/interPDB_analysis/Distance")

############################################################################################

# Read in df
df <- read.csv("Output/unique_chains_alignment.csv")

############################################################################################

# Inter-PDB RMSD

cat("Calculating Inter-PDB RMSD\n")

# Getting amyloid names from first chain only
amyloid_names <- unique(df$pdb_id)

#Initializing data frame to store mean_distance scores for each comparison
mean_distance_data <- data.frame(pdb_id = unique(df$pdb_id))


for (i in amyloid_names) {

  # Remove any positions where reference chain is not found (Unable to compare residues that don't exist in reference chain)
  filtered_df <- df %>%
                  group_by(residue_position) %>%
                    filter(any(pdb_id == i))
                
  # Creates new columns that equals the coordinates of the reference chain for each pdb
  filtered_df <- filtered_df %>%
    group_by(residue_position) %>%
    mutate(x_comp = x[pdb_id == i],
           y_comp = y[pdb_id == i],
           z_comp = z[pdb_id == i])
  
  # Finding distance between two points
  filtered_df <- filtered_df %>% 
                  mutate(distance = sqrt((x_comp - x)^2 + (y_comp - y)^2 + (z_comp - z)^2))
  
  
  # Plotting distance 
  
  ggplot(filtered_df, aes(x = residue_position, y = distance, colour = pdb_id)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1.2) +
    ggtitle(paste(i, " - Chain Comparison", sep = "")) +
    ylab("Distance") +
    xlab("Residue Position") +
    scale_x_continuous(breaks = seq(0,140,2)) +
    theme(panel.grid.major = element_line(colour = "grey30", linewidth = 0.5, 
                                          linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
          plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
          axis.text.y = element_text(size = 16, colour = "black"),
          axis.text.x = element_text(size = 16, colour = "black"),
          axis.title=element_text(size = 25, face = "bold"),
          axis.line.x = element_line(color="black", linewidth = 1.5),
          axis.line.y = element_line(color="black", linewidth = 1.5),
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 10, face = "bold")) 
  
  ggsave(plot = last_plot(), file = paste("Output/interPDB_analysis/Distance/", i, "_distance.png", sep = ""), 
         width = 20, height = 8)
  
  #################################################################################################
  #################################################################################################
  #################################################################################################
  #################################################################################################
  
  # Creating data frame for mean_distance heatmap
  
  # Finding the mean distance (i.e mean_distance) for each pdb
  comp_mean_distance <- filtered_df %>% group_by(pdb_id) %>% summarise_at(vars(distance), funs(mean))
  
  # Plotting comp_mean_distance 
  
  ggplot(comp_mean_distance, aes(x = reorder(pdb_id, distance), y = distance)) +
    geom_point(size = 4) +
    ggtitle(paste(i, " - Chain Comparison", sep = "")) +
    ylab("RMSD") +
    xlab("PDB") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(comp_mean_distance$distance) * 1.05)) +
    theme(panel.grid.major = element_line(colour = "grey30", linewidth = 0.5, 
                                          linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
          plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5), 
          axis.title = element_text(size = 18, face = "bold", colour = "black"),
          axis.text.x = element_text(size = 14, colour = "black",
                                     angle = 90, vjust = 0.4, hjust = 1),
          axis.text.y = element_text(size = 14, colour = "black"))
  
  ggsave(plot = last_plot(), file = paste0("Output/interPDB_analysis/RMSD/single_reference/", i, "_RMSD.png"), width = 14, height = 8)
  
  # Adding mean_distance to data frame
  mean_distance_data <- merge(mean_distance_data, comp_mean_distance, all = TRUE)
  
  # renaming last column so I know what comparison the similarities refer to
  names(mean_distance_data)[length(names(mean_distance_data))] <- i  
  
}

  
# Saving Mean Chain Distance Data 
#write.csv(mean_chain_distance_data, paste("Output/Figures/RMSD/",i, "_chain_comp_data.csv", sep = ""), row.names=FALSE)

# RMSD Heatmap Plotting

# RMSD = Root Mean Squared Distance

# Converting dataframe into long format for plotting
mean_distance_heatmap <- mean_distance_data %>% pivot_longer(!pdb_id, names_to = "comparison", values_to = "mean_distance")

# Heatmap plot
ggplot(mean_distance_heatmap, aes(x = comparison, y = pdb_id, fill= mean_distance)) + 
  geom_tile() +
  ggtitle("RMSD Heatmap") +
  ylab("PDB") +
  xlab("Reference PDB") +
  scale_fill_gradient(low="black", high="red", na.value = "white") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle = 90, 
                                   vjust = 0.4, hjust = 1),
        axis.title=element_text(size = 25, face = "bold"),
        axis.line.x = element_line(color="black", linewidth = 1.5),
        axis.line.y = element_line(color="black", linewidth = 1.5),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) +
  labs(fill = "RMSD")

# Saving the heatmap
ggsave(plot = last_plot(), file = paste("Output/interPDB_analysis/RMSD/RMSD_heatmap.png", sep = ""), 
       width = 14, height = 12)

# Saving data
write.csv(mean_distance_data, file = "Output/interPDB_analysis/RMSD/mean_distance_data.csv", row.names = FALSE)
write.csv(mean_distance_heatmap, file = "Output/interPDB_analysis/RMSD/mean_distance_heatmap.csv", row.names = FALSE)


# Cleaning Environment
#rm(list = ls())
invisible(gc())



