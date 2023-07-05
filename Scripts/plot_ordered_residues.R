# Installing and Loading Packages
library("readr")
library("dplyr")
library("stringr")
library("tidyr")
library("purrr")
library("ggplot2")

# Read in PDB metadata
data <- read_delim("Output/selected_pdbs_metadata.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

# Selecting only PDB and Residues columns
residues <- subset(data, select = c("PDB ID", "Residues Ordered"))

# Splitting Residues Ordered by "," as these indicate gaps within the chain
residues$`Residues Ordered` <- str_split(residues$`Residues Ordered`, ",")

# Assigning a unique group ID to each stretch of continuous residues in the pdb
# If for one PDB residues 1-10 and 20-60 are solved they will be given different groups
# This is so that the plot can draw lines to continuous regions and leave spaces where residues aren't ordered
residues <- residues %>%
              mutate(dash_count = str_count(`Residues Ordered`, "-"),
                     group = map2(dash_count, 1:dash_count, seq), 
                     group = str_split(group, ":")) %>%
              select(`PDB ID`, `Residues Ordered`, group)

# Unnesting the split character vectors into seperate rows and removing c("") symbols left over
residues <- residues %>% 
              unnest(cols = c(`Residues Ordered`, group)) %>% 
                mutate(`Residues Ordered` = str_replace_all(`Residues Ordered`, 
                                                             c("c" = "", "\"" = "", "\\(" = "", "\\)" = "")))

# Splitting each continuous residue group to give a min and max range
# e.g. 20-96 gets split into two rows; 20 and 96 
# On the plot a line will be draw connecting these two dots
residues$`Residues Ordered` <- str_split(residues$`Residues Ordered`, "-")

# Unnesting into seperate rows
df_unnested <- residues %>% unnest_longer(`Residues Ordered`)

# Setting `Residues Ordered` to be numeric
df_unnested$`Residues Ordered` <- as.numeric(df_unnested$`Residues Ordered`)

################################################################################

# Plotting
p <- ggplot(df_unnested, aes(x = `PDB ID`, y = `Residues Ordered`)) 
  
  if (any(grepl("synuclein", data$Protein))) {
    
    q <- p + # N-term
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 60, 
               alpha = 0.5, fill =  "darkblue", colour = "black") +
      # NAC
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 61, ymax = 95, 
               alpha = 0.5, fill =  "black", colour = "black") +
      # C-term
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 96, ymax = 140, 
               alpha = 0.5, fill =  "red", colour = "black") +
      # P1
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 36, ymax = 42, 
               alpha = 0.5, fill =  "orange", colour = "black") +
      # P2
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 45, ymax = 57, 
               alpha = 0.5, fill =  "green", colour = "black")
    
  } else if (any(grepl("Amyloid-", data$Protein))) {
    
    # Reference - https://doi.org/10.1021/jp210019h
    
    q <- p + # N-Term
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 16, 
               alpha = 0.5, fill =  "darkblue", colour = "black") +
      # Central Hydrophobic Core
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 17, ymax = 21, 
               alpha = 0.5, fill =  "black", colour = "black") +
      # Turn Region
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 24, ymax = 27, 
               alpha = 0.5, fill =  "red", colour = "black") +
      # Second Hydrophobic Region
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 28, ymax = 35, 
               alpha = 0.5, fill =  "orange", colour = "black") +
      # C-Term
      annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 36, ymax = 42, 
               alpha = 0.5, fill =  "green", colour = "black")
      
  } else { 
    q <- p 
  }
  
r <- q + geom_point(size = 3) +
  geom_line(aes(group = `PDB ID`), linewidth = 1, linetype = "dotted", colour = "black") +
  geom_line(aes(group = interaction(`PDB ID`, group)), linewidth = 2) +
  ggtitle("Ordered Residues") +
  ylab("Residue Position") +
  xlab("PDB") +
  scale_y_continuous(breaks = seq(0,140,5)) +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
        plot.title = element_text(size = 20, face = "bold", colour = "black", hjust = 0.5),
        axis.text.y = element_text(size = 12, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 12, colour = "black", face = "bold", angle = 45, hjust=1),
        axis.title=element_text(size = 25, face = "bold"),
        axis.line.x = element_line(color="black", linewidth = 1.5),
        axis.line.y = element_line(color="black", linewidth = 1.5),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 10, face = "bold")) 

plot(r)

# Saving Plot
ggsave(plot = last_plot(), file = "Output/Ordered_Residues.png", height = 6, width = 14)




