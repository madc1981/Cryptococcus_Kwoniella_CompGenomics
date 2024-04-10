#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: TE Classes Distribution in Cryptococcus and Kwoniella Species
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-05-10
Description: This script visualizes the distribution of Transposable Elements (TE) classes across Cryptococcus and
Kwoniella species, following a specific order for species and TE categories. It filters the input data 
to exclude specific TE categories, pivots the data for visualization, and generates a stacked bar plot to illustrate 
the TE class distribution for each species.

Requirements: pandas, matplotlib, seaborn
Input file: Cryptococcus_Kwoniella_EarlGrey_TEs_merged.txt
Output file: Cryptococcus_Kwoniella_TE_Classes_Distribution.pdf
---------------------------------------------------------------------------------------------------------------------
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
file_path = 'Cryptococcus_Kwoniella_EarlGrey_TEs_merged.txt'
data = pd.read_csv(file_path, sep='\t')

# Define organism order for the plot, then reverse the DataFrame rows before plotting
species_order = [
    "Cryptococcus_neoformans_H99",
    "Cryptococcus_deneoformans_JEC21",
    "Cryptococcus_bacillisporus_CA1280",
    "Cryptococcus_decagattii_7685027",
    "Cryptococcus_gatti_WM276",
    "Cryptococcus_tetragattii_IND107",
    "Cryptococcus_sp._MF34",
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
    "Kwoniella_shandongensis_CBS12478"
]

# Filtering and preparing the data
filtered_data = data[~data['tclassif'].isin(['Non-Repeat', 'Penelope', 'SINE'])]
stacked_data_filtered = filtered_data.pivot_table(index='Organism', columns='tclassif', values='percentage', fill_value=0)
stacked_data_filtered = stacked_data_filtered.reindex(species_order).iloc[::-1]

# Update the category name for simplicity
stacked_data_filtered.rename(columns={'Other (Simple Repeat, Microsatellite, RNA)': 'Others'}, inplace=True)

# Define the new order with the simplified category name
te_category_order = ['Unclassified', 'Rolling Circle', 'Others', 'LTR', 'LINE', 'DNA']

# Make sure to filter the DataFrame with the updated category names
stacked_data_filtered = stacked_data_filtered[te_category_order]

# Update values to reflect percentages
stacked_data_filtered *= 100

# Colors for TE categories
colors = ['#bcbec0', '#ff7f00', '#984ea3', '#4daf4a', '#377eb8', '#e41a1c']

# Plotting adjustments
plt.figure(figsize=(10, 14))  # Adjust figure size as needed
bar_width = 0.8  # Bar width

# Stacked bar plot
ax = stacked_data_filtered.plot(kind='barh', stacked=True, color=colors, edgecolor='none', width=bar_width)

# Starts x-axis slightly before zero to give space
ax.set_xlim(left=-0.3)

plt.title('TE Classes Distribution by Species', fontsize=10)
plt.xlabel('Percentage', fontsize=10)
plt.ylabel('Species', fontsize=10)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

# Adding dashed lines to the x-axis scale bar covering the full plot
ax.set_axisbelow(True)  # Ensures grid lines are plotted below the bars
ax.xaxis.grid(color='gray', linestyle='dashed', linewidth=0.5)

# Adjusting subplot parameters
plt.subplots_adjust(top=0.928, bottom=0.113, left=0.417, right=0.812, hspace=0.2, wspace=0.2)

# Update legend placement and title
handles, labels = ax.get_legend_handles_labels()
labels = [label if label != 'Other (Simple Repeat, Microsatellite, RNA)' else 'Others' for label in labels]

# Optimal legend placement inside the plot area
ax.legend(handles, labels, title='Transposon class', loc='lower right', fontsize=8)

plt.rcParams['pdf.fonttype'] = 42
plt.savefig('Cryptococcus_Kwoniella_TE_Classes_Distribution.pdf', bbox_inches='tight')
plt.show()