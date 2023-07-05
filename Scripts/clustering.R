cat("Running RMSD_analysis.R\n")

# Installing and Loading Packages
library("dplyr")
library("stringr")
library("tidyr")
library("readr")
library("ggplot2")
library("dendextend")

############################################################################################

# Creating Directories
if(!dir.exists("Output/Clustering")) dir.create("Output/Clustering")

############################################################################################

# Read in data
chain_data <- read_csv("Output/single_chain_analysis/unique_chain_data.csv")

# Creating pdb_id so if PDB has multiple distinct chains they are saved as PDB_chain
chain_data$pdb_id <- ""

for (i in 1:nrow(chain_data)) {
  
  if (chain_data$chain[i] == "A") {
    
    chain_data$pdb_id[i] <- as.character(chain_data$pdb[i])
    
  } else {
    
    chain_data$pdb_id[i] <- paste0(as.character(chain_data$pdb[i]), "_", as.character(chain_data$chain[i]))
    
  }
  
}

# Importing RMSD Data
RMSD_data <- read_csv("Output/interPDB_analysis/RMSD/mean_distance_data.csv")

# Selecting only 1 RMSD Comparison for clustering so it doesn't dominate comparison metrics

# Getting .pdb names for RMSD Calculations
pdb_names <- list.files(pattern = "\\.pdb$")
pdb_names <- gsub("\\.pdb$", "", pdb_names)

if (length(pdb_names) > 0) {

  # Selecting RMSD values from first pdb comparison
  RMSD_data <- subset(RMSD_data, select = c("pdb_id", pdb_names[1]))
  
} else {
  
  # If no pdb files, select first RMSD column
  RMSD_data <- RMSD_data[,1:2]
  
}


df <- merge(chain_data, RMSD_data, by = "pdb_id") 

# Find NA values to remove from heatmap data later
RMSD_NA <- subset(df, is.na(df[, ncol(df)]), select = pdb_id)

# Remove rows with NA values
df <- na.omit(df)

############################################################################################

# Removing reference columns
cluster_data <- subset(df, select = setdiff(names(df), c("pdb", "chain", "pdb_id")))

# Scaling Cluster Data (makes mean of all columns = 0 and SD = 0-1)
scaled_cluster_data <- as.data.frame(scale(cluster_data))
summary(scaled_cluster_data)

# Calculating Distance
cluster_dist <- dist(scaled_cluster_data, method = "euclidean")

# Clustering
hc_average <- hclust(cluster_dist, method = "average")

# Making labels
pdb_id_labels <- df$pdb_id[hc_average$order]

# Cutting Dendrogram to Create Clusters

# Setting cut_height to the mean distance
cut_height <- mean(cluster_dist)

# Extracting pdb_id from each cluster group
clusters <- cutree(hc_average, h = cut_height)
grouped_pdb_ids <- split(df$pdb_id, clusters)

# Create df of each pdb_id and its associated group
grouping_result <- data.frame(pdb_id = unlist(grouped_pdb_ids), group = rep(1:length(grouped_pdb_ids), sapply(grouped_pdb_ids, length)))

# Splitting pdb_id into PDB and chain
grouping_result$pdb <- ""
grouping_result$chain <- "A"

for (i in 1:nrow(grouping_result)) {

  grouping_result$pdb[i] <- as.character(grouping_result$pdb_id[i])
  
  if (str_sub(grouping_result$pdb_id[i], start = -2, end = -2) == "_") {
    
    grouping_result$pdb[i] <- strsplit(as.character(grouping_result$pdb_id[i]), "_")[[1]][1]
    
    grouping_result$chain[i] <- str_sub(as.character(grouping_result$pdb_id[i]), start = -1, end = -1)
    
  } 
    
}

# Reordering columns
grouping_result <- grouping_result[, c("pdb", "chain", "pdb_id", "group")]

# Saving grouping result
write.csv(grouping_result, file = "Output/Clustering/grouped_pdb_ids.csv", row.names = FALSE)

# Changing line colours instead of adding borders
avg_dend_obj <- as.dendrogram(hc_average)
avg_col_dend <- color_branches(avg_dend_obj, h = cut_height)
labels(avg_col_dend) <- pdb_id_labels
plot(avg_col_dend, ylab = "Euclidean Distance") # horiz = TRUE

# Save the cluster plot as a PNG file
png("Output/Clustering/cluster_dendrogram.png", width = 4000, height = 2500, res = 300)
plot(avg_col_dend, ylab = "Euclidean Distance")
dev.off()

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

# Creating an RMSD heatmap for all chains and ordering it by cluster

# importing heatmap plot data
mean_distance_heatmap <- read_csv("Output/interPDB_analysis/RMSD/mean_distance_heatmap.csv")

# Removing PDBs with NAs
mean_distance_heatmap <- subset(mean_distance_heatmap, !(pdb_id %in% RMSD_NA$pdb_id | comparison %in% RMSD_NA$pdb_id))


# Ordering by cluster order
mean_distance_heatmap$pdb_id <- factor(mean_distance_heatmap$pdb_id, levels = pdb_id_labels)
mean_distance_heatmap$comparison <- factor(mean_distance_heatmap$comparison, levels = pdb_id_labels)

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
ggsave(plot = last_plot(), file = paste("Output/Clustering/RMSD_Heatmap.png", sep = ""), 
       width = 14, height = 12)




