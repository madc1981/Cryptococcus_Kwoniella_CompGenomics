# Script Metadata
#------------------------------------------------------------------------------------------
# Title: Ridgeline Plot for GC Content in Cryptococcus and Kwoniella Species
# Author: Marco A. Coelho @ Heitman lab
# Date: 2023-04-06
# Description: This script generates a ridgeline plot visualizing GC content distribution
#              across Cryptococcus and Kwoniella species. It leverages ggplot2
#              and ggridges for visualization, and dplyr for data manipulation.
# Input file:  "Cryptococcus_Kwoniella_GCcontent.headers.txt"
#------------------------------------------------------------------------------------------

# Load necessary libraries
library(ggplot2) # For creating plots
library(dplyr) # For data manipulation
library(ggridges) # For creating ridgeline plots

# Data Preparation
#------------------------------------------------------------------------------------------
# Note: It's generally recommended to use relative paths or set the working directory outside the script for better reproducibility
# Example: setwd("./data")

# Read the data
df <- read.table("Cryptococcus_Kwoniella_GCcontent_1kb_window.txt", header = TRUE)

# Calculate the mean GC content for each species
mu <- df %>%
  group_by(Species) %>%
  summarise(grp.mean = mean(GCperc))

# Define the order of species to display in the plot
species_order <- rev(c("Cryptococcus_neoformans_H99",
                       "Cryptococcus_deneoformans_JEC21",
                       "Cryptococcus_bacillisporus_CA1280",
                       "Cryptococcus_decagattii_7685027",
                       "Cryptococcus_gatti_WM276",
                       "Cryptococcus_tetragattii_IND107",
                       "Cryptococcus_gatti_MF34",
                       "Cryptococcus_deuterogattii_R265",
                       "Cryptococcus_depauperatus_CBS7841",
                       "Cryptococcus_floricola_DSM27421",
                       "Cryptococcus_wingfieldii_CBS7118",
                       "Cryptococcus_amylolentus_CBS6039",
                       "Cryptococcus_sp._OR849",
                       "Cryptococcus_sp._OR918",
                       "Kwoniella_europaea_PYCC6329",
                       "Kwoniella_sp._B9012",
                       "Kwoniella_botswanensis_CBS12716",
                       "Kwoniella_mangrovensis_CBS8507",
                       "Kwoniella_bestiolae_CBS10118",
                       "Kwoniella_dejecticola_CBS10117",
                       "Kwoniella_pini_CBS10737",
                       "Kwoniella_dendrophila_CBS6074",
                       "Kwoniella_shivajii_CBS11374",
                       "Kwoniella_heveanensis_CBS569",
                       "Kwoniella_sp._CBS6097",
                       "Kwoniella_sp._CBS9459",
                       "Kwoniella_sp._DSM27419",
                       "Kwoniella_newhampshirensis_CBS13917",
                       "Kwoniella_shandongensis_CBS12478"))

# Set the factor levels for Species based on the desired order for consistent plotting
df$Species <- factor(df$Species, levels = species_order)

# Visualization
#------------------------------------------------------------------------------------------
# Color mappings for both fill and line colors reflecting group affiliations
fill_color_mapping <- c(rep("#bae4f6", 15), rep("#f9bfb2", 14)) # Fill colors for Kwoniella and Cryptococcus, respectively
line_color_mapping <- c(rep("#6f98a9", 15), rep("#ac766a", 14)) # Line colors for Kwoniella and Cryptococcus, respectively
names(fill_color_mapping) <- species_order
names(line_color_mapping) <- species_order

# Generate and customize the Ridgeline plot
my_plot <- ggplot(df, aes(x = GCperc, y = Species, fill = Species, color = Species)) + 
  geom_density_ridges(rel_min_height = 0.01) +
  scale_fill_manual(values = fill_color_mapping) +
  scale_color_manual(values = line_color_mapping) +
  xlim(15, 75) + # Restricting x-axis range to focus on relevant GC content percentages
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 14) +
  stat_density_ridges(
    geom = "density_ridges_gradient", quantile_lines = TRUE, calc_ecdf = TRUE, 
    alpha = 0.1, quantiles = 2,
    position = "identity"
  ) +
  theme(legend.position = "none")

# Display the plot
print(my_plot)

# Saving the Plot
#------------------------------------------------------------------------------------------
# Save the generated plot as a PDF, specifying desired dimensions and filename
ggsave("Cryptococcus_Kwoniella_GCcontent_ridgeline_plot.pdf", plot = my_plot, width = 6, height = 12)