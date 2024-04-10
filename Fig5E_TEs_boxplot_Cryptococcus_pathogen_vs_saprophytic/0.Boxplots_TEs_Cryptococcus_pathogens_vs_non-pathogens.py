#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: Comparative Analysis of Transposable Element Content in pathogen vs. saprophytic Cryptococcus Species
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-05-09
Description: This script analyzes the differences in transposable element (TE) content between
pathogen and saprophytic Cryptococcus species. It filters input data to exclude certain
TE categories, calculates total TE percentages for each species, classifies species based on
Group, and performs statistical tests to compare TE content. The results are visualized
in a boxplot, annotated with statistical test results, and saved along with a statistical summary.

Requirements: pandas, matplotlib, seaborn, scipy
Input file: Cryptococcus_and_Kwoniella_EarlGrey_TEs_merged.txt
Output files: TEs_Cryptococcus_pathogens_vs_Cryptococcus_non-pathogens.pdf,
              TEs_Cryptococcus_pathogens_vs_Cryptococcus_non-pathogens_Stats.csv
              TEs_Cryptococcus_pathogens_vs_saprophytic_Table.csv
---------------------------------------------------------------------------------------------------------------------
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# Load the data
# Load the data
file_path = 'Cryptococcus_Kwoniella_EarlGrey_TEs_merged.txt'
data = pd.read_csv(file_path, sep='\t')

# Define pathogen and saprophytic Cryptococcus species
pathogen_species = [
    'Cryptococcus_neoformans_H99',
    'Cryptococcus_deneoformans_JEC21',
    'Cryptococcus_bacillisporus_CA1280',
    'Cryptococcus_decagattii_7685027',
    'Cryptococcus_gatti_WM276',
    'Cryptococcus_tetragattii_IND107',
    'Cryptococcus_sp._MF34',
    'Cryptococcus_deuterogattii_R265'
]

saprophytic_species = [
    'Cryptococcus_depauperatus_CBS7841',
    'Cryptococcus_floricola_DSM27421',
    'Cryptococcus_wingfieldii_CBS7118',
    'Cryptococcus_amylolentus_CBS6039',
    'Cryptococcus_sp._OR849',
    'Cryptococcus_sp._OR918',
]

# Filter out TE categories
boxplot_data = data[~data['tclassif'].isin(['Non-Repeat', 'Penelope', 'SINE'])]

# Calculate total percentage of transposons for each organism
total_transposons = boxplot_data.groupby('Organism')['percentage'].sum().reset_index()

# Classify each organism as pathogen or saprophytic
total_transposons['Group'] = total_transposons['Organism'].apply(
    lambda x: 'pathogen' if x in pathogen_species else ('saprophytic' if x in saprophytic_species else 'Other'))

# Filter out organisms not in the specified lists
total_transposons = total_transposons[total_transposons['Group'] != 'Other']

# Adjusting the 'percentage' values to reflect true percentages
total_transposons['percentage'] *= 100

# Combine pathogen and saprophytic species into one ordered list
combined_species_order = pathogen_species + saprophytic_species

# Create a sorting key based on the combined order
total_transposons['Sort_Key'] = total_transposons['Organism'].apply(lambda x: combined_species_order.index(x))

# Sort the DataFrame by the custom order
total_transposons_sorted = total_transposons.sort_values('Sort_Key').drop('Sort_Key', axis=1)

# Save the sorted raw data used for the boxplot
total_transposons_sorted.to_csv('TEs_Cryptococcus_pathogens_vs_saprophytic_Table.csv', index=False)


# Boxplot with mean value and data points
plt.figure(figsize=(4, 8))
sns.boxplot(x='Group', y='percentage', data=total_transposons, 
            order=['pathogen', 'saprophytic'],  # Specifying the order here
            showmeans=True, meanline=True, 
            palette={'pathogen': "#fad2a4", 'saprophytic': "#9ad3b5"})
sns.stripplot(x='Group', y='percentage', data=total_transposons, 
              order=['pathogen', 'saprophytic'],  # Specifying the order here
              color='grey', edgecolor='black', size=10, jitter=True)

# Statistical Tests
pathogen_transposons = total_transposons[total_transposons['Group'] == 'pathogen']['percentage']
saprophytic_transposons = total_transposons[total_transposons['Group'] == 'saprophytic']['percentage']

# Shapiro-Wilk Test for normality and Mann-Whitney U Test
p_pathogen = stats.shapiro(pathogen_transposons).pvalue
p_non_pathogen = stats.shapiro(saprophytic_transposons).pvalue
u_stat, u_p_value = stats.mannwhitneyu(pathogen_transposons, saprophytic_transposons)

# Plot annotations with word-wrapped title
plt.title('Percentage of Transposons in\npathogens vs. saprophytic Cryptococcus', pad=10)
plt.xlabel('Group')
plt.ylabel('Percentage of TEs (%)')  # Updated label to reflect percentage

# Adjusting subplot parameters and saving the plot as PDF
plt.subplots_adjust(top=0.92, bottom=0.06, left=0.157, right=0.969, hspace=0.2, wspace=0.2)
plt.rcParams['pdf.fonttype'] = 42
plt.savefig('TEs_Cryptococcus_pathogens_vs_saprophytic.pdf')
plt.show()

# Find the species with the min and max percentages for both groups
min_pathogen_species = total_transposons[total_transposons['Group'] == 'pathogen'].sort_values('percentage').iloc[0]['Organism']
max_pathogen_species = total_transposons[total_transposons['Group'] == 'pathogen'].sort_values('percentage').iloc[-1]['Organism']
min_non_pathogen_species = total_transposons[total_transposons['Group'] == 'saprophytic'].sort_values('percentage').iloc[0]['Organism']
max_non_pathogen_species = total_transposons[total_transposons['Group'] == 'saprophytic'].sort_values('percentage').iloc[-1]['Organism']

# Preparing statistical results for CSV
stats_results = {
    'Statistic': ['Mean pathogen', 'Mean saprophytic', 'Median pathogen', 'Median saprophytic',
                  'Min pathogen', 'Min saprophytic', 'Max pathogen', 'Max saprophytic',
                  'Shapiro-Wilk pathogen', 'Shapiro-Wilk saprophytic', 'Mann-Whitney U'],
    'Value': [pathogen_transposons.mean(), saprophytic_transposons.mean(), 
              pathogen_transposons.median(), saprophytic_transposons.median(),
              pathogen_transposons.min(), saprophytic_transposons.min(), 
              pathogen_transposons.max(), saprophytic_transposons.max(),
              None, None, u_stat],
    'p-value': [None, None, None, None, None, None, None, None,
                p_pathogen, p_non_pathogen, u_p_value],
    'Species': [None, None, None, None, min_pathogen_species, min_non_pathogen_species, 
                max_pathogen_species, max_non_pathogen_species, None, None, None]
}

stats_df = pd.DataFrame(stats_results)
stats_df.to_csv('TEs_Cryptococcus_pathogens_vs_saprophytic_Stats.csv', index=False)

print('Analysis complete. Files saved:')
print('TEs_Cryptococcus_pathogens_vs_saprophytic.pdf')
print('TEs_Cryptococcus_pathogens_vs_saprophytic_Stats.csv')
print('TEs_Cryptococcus_pathogens_vs_saprophytic_Table.csv')
