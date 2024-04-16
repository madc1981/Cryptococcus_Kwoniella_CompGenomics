#!/usr/bin/env Rscript

# Metadata
# ---------------------------------------------------------------------------------------------------------------------
# Title: UpSet plot for Orthogroups (OGs) categorization
# Author: Marco A. Coelho @ Heitman lab, Duke University
# Date: 2023-04-15
# Description: This script generates an UpSet plot to visualize the intersection of orthogroups across multiple species
#              in a genomics dataset. It allows exclusion of specific species, defines specific groupings for visualization,
#              and highlights particular species groups in the plot.
# Input files: Orthogroups_renamed_sorted.tsv (renamed and sorted OrthoFinder file using the python script
#              0.rename_and_sort_OG_tables.py)
# Requirements: UpSetR
# ---------------------------------------------------------------------------------------------------------------------

# Install and load necessary packages if not already installed
if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)

# Read in the data
data <- read.delim("Orthogroups_renamed_sorted.tsv", header = TRUE, stringsAsFactors = FALSE)

# Convert data to a binary matrix
binary_data <- data
for (i in 2:ncol(data)) {
  binary_data[, i] <- ifelse(data[, i] != "", TRUE, FALSE)
}
# Exclude the Orthogroup column to only have the binary matrix
binary_data$Orthogroup <- NULL

# Exclude specified species from analysis (outgroups)
exclude_species <- c("Sait_podz_DSM27192", "Trem_mese_ATCC28783", "Bull_alba_JCM2954")
cols_to_exclude <- colnames(binary_data)[colnames(binary_data) %in% exclude_species]
if (length(cols_to_exclude) > 0) {
  binary_data <- binary_data[, !(colnames(binary_data) %in% cols_to_exclude)]
}

# Convert binary matrix to long format
long_data <- data.frame(
  Orthogroup = rep(rownames(binary_data), times = ncol(binary_data)),
  Species = factor(rep(colnames(binary_data), each = nrow(binary_data))),
  Presence = as.logical(as.vector(as.matrix(binary_data)))
)

# Filter for presence of orthogroups
filtered_data <- long_data[long_data$Presence == TRUE, ]

# Convert filtered long data to list format for UpSet plot
list_data <- split(filtered_data$Orthogroup, filtered_data$Species)

# Define species orders and group specifications for visualization
species_order <- colnames(rev(binary_data))
crypto_pathogens <- c("Cryp_neof_H99", "Cryp_dene_JEC21", "Cryp_baci_CA1280", "Cryp_deca_7685027", "Cryp_gatt_WM276", "Cryp_tetr_IND107", "Cryp_sp._MF34", "Cryp_deut_R265")
removing_Cdep <- c("Cryp_depa_CBS7841")
only_OR918 <- c("Cryp_sp._OR918")
Kwoniella_specific <- c("Kwon_euro_PYCC6329", "Kwon_sp._B9012", "Kwon_bots_CBS12716", "Kwon_mang_CBS8507", "Kwon_best_CBS10118", "Kwon_deje_CBS10117", "Kwon_pini_CBS10737", "Kwon_dend_CBS6074", "Kwon_shiv_CBS11374", "Kwon_heve_CBS569", "Kwon_sp._CBS6097", "Kwon_sp._CBS9459", "Kwon_sp._DSM27419", "Kwon_newh_CBS13917", "Kwon_shan_CBS12478")
only_Crypto <- c("Cryp_neof_H99", "Cryp_dene_JEC21", "Cryp_baci_CA1280", "Cryp_deca_7685027", "Cryp_gatt_WM276", "Cryp_tetr_IND107", "Cryp_sp._MF34", "Cryp_deut_R265", "Cryp_depa_CBS7841", "Cryp_flor_DSM27421", "Cryp_wing_CBS7118", "Cryp_amyl_CBS6039", "Cryp_sp._OR849", "Cryp_sp._OR918")
non_patho <- c("Cryp_neof_H99", "Cryp_dene_JEC21", "Cryp_baci_CA1280", "Cryp_deca_7685027", "Cryp_gatt_WM276", "Cryp_tetr_IND107", "Cryp_sp._MF34", "Cryp_deut_R265")

# Define queries for the UpSet plot
queries <- list(
  list(
    query = intersects,
    params = list(crypto_pathogens),
    active = TRUE,
    color = "#e17f35",
    query.name = "Cryptococcus pathogens"
  ),
  list(
    query = intersects,
    params = list(species_order[!species_order %in% removing_Cdep]),
    active = TRUE,
    color = "#88d392",
    query.name = "Losses in C. depauperatus"
  ),
  list(
    query = intersects,
    params = list(species_order[species_order == only_OR918]),
    active = TRUE,
    color = "#9b8579",
    query.name = "Cryptococcus sp. OR918 specific"
  ),
  list(
    query = intersects,
    params = list(Kwoniella_specific),
    active = TRUE,
    color = "#6da7bf",
    query.name = "Kwoniella specific"
  ),
  list(
    query = intersects,
    params = list(only_Crypto),
    active = TRUE,
    color = "#d87881",
    query.name = "Cryptococcus specific"
  ),
  list(
    query = intersects,
    params = list(species_order[!species_order %in% non_patho]),
    active = TRUE,
    color = "#918dc2",
    query.name = "Non-pathogens specific"
  )
)

# Generate the UpSet plot
upset(
  fromList(list_data),
  sets = species_order,
  query.legend = "bottom",
  queries = queries,
  keep.order = TRUE,
  mb.ratio = c(0.50, 0.50),
  order.by = "freq",
  nintersects = 45,
  nsets = length(list_data),
  number.angles = 30,  # Adjusts the angle for set names
  point.size = 3.5,
  line.size = 1.0,
  mainbar.y.label = "Intersected OGs",
  sets.x.label = "Number of OGs",
  # Changes the text size for different components
  # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  text.scale = c(2,1.5,2,1.5,1.5,1)
)

