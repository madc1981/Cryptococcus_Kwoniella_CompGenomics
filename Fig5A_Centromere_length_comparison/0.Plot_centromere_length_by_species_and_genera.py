
#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: Comparative analysis of centromere lengths in Cryptococcus and Kwoniella
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-05-10
Description: This script generates plots centromere lengths across Cryptococcus and Kwoniella species. Each centromere
             is represented by a dot, with median (dashed) and mean (solid) lines indicating typical values for each
             species. An inset provides a statistical comparison of centromere lengths between the two genera.
             Shapiro-Wilk tests assess normality, and appropriate statistical tests (Mann-Whitney U or t-test) compare
             the groups, with results exported to a CSV file and annotated on the inset plot.
Requirements: pandas, matplotlib, numpy, scipy, mpl_toolkits.axes_grid1

Input file: Cryptococcus_Kwoniella_Cen_length.txt
Output file: Cen_Length_Distribution_with_Inset.pdf
             Centromere_length_stats.csv
             Statistical_tests_results.csv
---------------------------------------------------------------------------------------------------------------------
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import shapiro, mannwhitneyu, ttest_ind
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Load and prepare the data
file_path = 'Cryptococcus_Kwoniella_Cen_length.txt'
data = pd.read_csv(file_path, sep='\t')
data['cen_length_kb'] = data['cen_length_bp'] / 1000
data['Group'] = data['Species'].apply(lambda x: 'Cryptococcus' if 'Cryptococcus' in x else 'Kwoniella')

# Correct species order as initially provided
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

# Split data into Cryptococcus and Kwoniella
cryptococcus_data = data[data['Group'] == 'Cryptococcus']['cen_length_kb']
kwoniella_data = data[data['Group'] == 'Kwoniella']['cen_length_kb']

# Perform Shapiro-Wilk test of normality
shapiro_cryptococcus = shapiro(cryptococcus_data)
shapiro_kwoniella = shapiro(kwoniella_data)

# Selecting statistical test based on normality test results
if shapiro_cryptococcus.pvalue < 0.05 or shapiro_kwoniella.pvalue < 0.05:
    # Non-normal distributions
    test_stat, p_value = mannwhitneyu(cryptococcus_data, kwoniella_data)
    test_name = "Mann-Whitney U"
else:
    # Normal distributions
    test_stat, p_value = ttest_ind(cryptococcus_data, kwoniella_data)
    test_name = "Independent Samples t-test"

# Format p-value for plot annotation
if p_value < 0.0001:
    p_value_str = "P < 0.0001"
else:
    p_value_str = f"P = {p_value:.4f}"

# Save the test results, including Shapiro-Wilk test results, to a CSV file
test_results_df = pd.DataFrame({
    'Test Name': ["Shapiro-Wilk Cryptococcus", "Shapiro-Wilk Kwoniella", test_name],
    'Statistic': [shapiro_cryptococcus.statistic, shapiro_kwoniella.statistic, test_stat],
    'p-value': [shapiro_cryptococcus.pvalue, shapiro_kwoniella.pvalue, p_value]
})
test_results_csv_path = 'Statistical_tests_results.csv'
test_results_df.to_csv(test_results_csv_path, index=False)

# Print the path to the CSV file with test results
print(f'Test results saved to {test_results_csv_path}')

# Initialize a list to store median and mean values for each species
centromere_stats = []

for species in species_order:
    species_data = data[data['Species'] == species]['cen_length_kb']
    median_val = species_data.median()
    mean_val = species_data.mean()
    centromere_stats.append({'Species': species, 'Median length (kb)': median_val, 'Mean length (kb)': mean_val})

# Calculate group statistics for Cryptococcus and Kwoniella
cryptococcus_group_stats = {
    'Species': 'Cryptococcus Group',
    'Median Length (kb)': cryptococcus_data.median(),
    'Mean Length (kb)': cryptococcus_data.mean()
}

kwoniella_group_stats = {
    'Species': 'Kwoniella Group',
    'Median Length (kb)': kwoniella_data.median(),
    'Mean Length (kb)': kwoniella_data.mean()
}

# Append group statistics to the list
centromere_stats.append(cryptococcus_group_stats)
centromere_stats.append(kwoniella_group_stats)

# Convert the list to a DataFrame and export to CSV as before
centromere_stats_df = pd.DataFrame(centromere_stats)

# Export the DataFrame to a CSV file
centromere_stats_csv_path = 'Centromere_Length_Stats.csv'
centromere_stats_df.to_csv(centromere_stats_csv_path, index=False)

# Print the path to the CSV file with median and mean centromere lengths
print(f'Centromere length statistics saved to {centromere_stats_csv_path}')
# Main plot
plt.figure(figsize=(8, 10))
ax_main = plt.gca()

for i, species in enumerate(species_order):
    species_data = data[data['Species'] == species]
    color = '#f7a38f' if 'Cryptococcus' in species else '#9bd8f2'
    
    # Plot each centromere as a dot
    ax_main.plot(species_data['cen_length_kb'], np.repeat(i, len(species_data)), 'o', color=color, markersize=10, alpha=0.6)
    
    # Plot median and mean as lines
    median_val = species_data['cen_length_kb'].median()
    mean_val = species_data['cen_length_kb'].mean()
    plt.plot([median_val]*2, [i-0.4, i+0.4], '--', color='black', linewidth=1.5)
    plt.plot([mean_val]*2, [i-0.4, i+0.4], '-', color='black', linewidth=1.5)

# Mean and median lines are plotted within the loop for individual groups

ax_main.set_yticks(range(len(species_order)))
ax_main.set_yticklabels(species_order, fontsize=8)
ax_main.invert_yaxis()  # Reverses the y-axis
ax_main.set_xlabel('Centromere length (kb)', fontsize=12)
ax_main.set_xlim(-3, 143)  # Extend x-axis and provide space at start
ax_main.tick_params(axis='x', labelsize=12)  # Adjust x-axis tick label fontsize
ax_main.tick_params(axis='y', labelsize=10)  # Adjust y-axis tick label fontsize

# Inset plot for comparison
# Creating inset within the main plot
ax_inset = inset_axes(ax_main, width="30%", height="30%", loc='lower right', borderpad=3)

# Plot individual points for Cryptococcus and Kwoniella
cryptococcus_points = data[data['Group'] == 'Cryptococcus']['cen_length_kb']
kwoniella_points = data[data['Group'] == 'Kwoniella']['cen_length_kb']

# Calculate mean and median values for Cryptococcus and Kwoniella
crypt_mean = cryptococcus_data.mean()
crypt_median = cryptococcus_data.median()
kwon_mean = kwoniella_data.mean()
kwon_median = kwoniella_data.median()

# Generate jittered x positions for Cryptococcus and Kwoniella points
cryptococcus_jitter = np.random.normal(loc=0, scale=0.1, size=len(cryptococcus_points))
kwoniella_jitter = np.random.normal(loc=1, scale=0.1, size=len(kwoniella_points))

# Plot points with jitter
ax_inset.plot(cryptococcus_jitter, cryptococcus_points, 'o', color='#f7a38f', markersize=5, label='Cryptococcus', alpha=0.6)
ax_inset.plot(kwoniella_jitter, kwoniella_points, 'o', color='#9bd8f2', markersize=5, label='Kwoniella', alpha=0.6)
ax_inset.set_xlim(-0.5, 1.5)  # Extend x-axis and provide space at start
ax_inset.set_ylim(-3, 143)  # Extend x-axis and provide space at start

# Mean and median lines for Cryptococcus and Kwoniella in the inset plot
ax_inset.axhline(crypt_median, color='black', linestyle='--', linewidth=1, xmin=0.1, xmax=0.4, label='Crypt. Median')
ax_inset.axhline(crypt_mean, color='black', linestyle='-', linewidth=1, xmin=0.1, xmax=0.4, label='Crypt. Mean')
ax_inset.axhline(kwon_median, color='black', linestyle='--', linewidth=1, xmin=0.6, xmax=0.9, label='Kwon. Median')
ax_inset.axhline(kwon_mean, color='black', linestyle='-', linewidth=1, xmin=0.6, xmax=0.9, label='Kwon. Mean')

# Adjust the x-ticks after adding the lines
ax_inset.set_xticks([0, 1])
ax_inset.set_xticklabels(['Crypt.', 'Kwon.'])
ax_inset.set_ylabel('Centromere length (kb)', fontsize=12)
ax_inset.tick_params(axis='x', labelsize=12)  # Adjust x-axis tick label fontsize for the inset plot
ax_inset.tick_params(axis='y', labelsize=12)  # Adjust y-axis tick label fontsize for the inset plot
plt.setp(ax_inset.get_xticklabels(), ha="center", rotation_mode="anchor")

# Adding p-value annotation to the inset plot
p_value_annotation = f'{test_name}\n {p_value_str}'
ax_inset.annotate(p_value_annotation, xy=(0.5, 1.05), xycoords='axes fraction', ha='center',
                  fontsize=10)

# Custom lines for the legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='black', lw=1.5),
                Line2D([0], [0], color='black', lw=1.5, linestyle='--')]

# Create the legend and add it to ax_main
ax_main.legend(custom_lines, ['Mean', 'Median'], loc='upper right', fontsize=10)

# Adjust subplot and save the plot
plt.tight_layout()
plt.rcParams['pdf.fonttype'] = 42
plt.savefig('Cen_Length_Distribution_with_Inset.pdf')
plt.show()
