#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------
Title: Comparative LTR TE distribution in Centromeric (CEN) vs Non-Centromeric (non-CEN) Regions
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-05-10
Description: This script visualizes the distribution of Long Terminal Repeat (LTR) transposable elements 
             in centromeric and non-centromeric regions across Cryptococcus and Kwoniella species. 
Requirements: pandas, matplotlib

Input file: Cryptococcus_Kwoniella_cen_and_non-cen_LTRs.txt
Output file: LTR_TE_Distribution_Centromeric_vs_NonCentromeric.pdf
---------------------------------------------------------------------------------------------------------------
"""

import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = 'Cryptococcus_Kwoniella_cen_and_non-cen_LTRs.txt'
data = pd.read_csv(file_path, sep='\t')

# Reverse the order of species for y-axis display
data = data.iloc[::-1]

# Plotting adjustments
plt.figure(figsize=(10, 14))
bar_width = 0.8  # Bar width

# Use the correct percentage columns for plotting
ax = data.set_index('Organism')[['normalized_TE_percent_cen', 'normalized_TE_percent_non_cen']].plot(
    kind='barh', stacked=True, color=['#af9d94', '#dad3c3'], edgecolor='none', width=bar_width)

# Starts x-axis slightly before zero to give space
ax.set_xlim(left=-3)

plt.title('LTRs distribution in \nCentromeric vs Non-Centromeric Regions', fontsize=10)
plt.xlabel('Percentage (%)', fontsize=10)
plt.ylabel('Species', fontsize=10)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.legend(['CEN', 'Non-CEN'], loc='lower right', fontsize=8)

# Adding dashed lines to the x-axis scale bar covering the full plot
ax.set_axisbelow(True)  # Ensures grid lines are plotted below the bars
ax.xaxis.grid(color='gray', linestyle='dashed', linewidth=0.5)

# Adjusting subplot parameters
plt.subplots_adjust(top=0.928, bottom=0.113, left=0.417, right=0.812, hspace=0.2, wspace=0.2)

plt.rcParams['pdf.fonttype'] = 42
plt.savefig('LTR_TE_distribution_CEN_vs_NonCEN_Cryptococcus_Kwoniella.pdf', bbox_inches='tight')
plt.show()
