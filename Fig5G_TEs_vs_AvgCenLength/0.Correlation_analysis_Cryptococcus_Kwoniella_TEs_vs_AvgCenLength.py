#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------
Title: Correlation Analysis of Avg. CEN length and Transposable Element (TE) Percentage
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-05-08
Description: This script performs a correlation analysis between Avg. CEN length and
             the percentage of transposable elements (TEs) in Cryptococcus and 
             Kwoniella species. It filters the input data to exclude certain TE 
             categories, calculates total TE percentages, assesses the data for 
             normal distribution, and decides on using Pearson or Spearman 
             correlation. The results are plotted, annotated, and saved alongside 
             the statistical and merged dataset outputs.
Requirements: pandas, matplotlib, seaborn, scipy, numpy
Input files: Cryptococcus_Kwoniella_EarlGrey_TEs_merged.txt
             Cryptococcus_Kwoniella_cen_length.txt
---------------------------------------------------------------------------------------
"""

# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

# Load the data
file_path = 'Cryptococcus_Kwoniella_EarlGrey_TEs_merged.txt'
data = pd.read_csv(file_path, sep='\t')

# Data Preprocessing
#--------------------------------------------------------------------------------
# Filter out specified TE categories
data_filtered = data[~data['tclassif'].isin(['Non-Repeat', 'Penelope', 'SINE'])]

# Calculate total TE percentage for each species
total_percentage = data_filtered.groupby('Organism')['percentage'].sum().reset_index()

# Load the centromere length data
file_path_centromere = 'Cryptococcus_Kwoniella_cen_length.txt'
centromere_data = pd.read_csv(file_path_centromere, sep='\t')

# Merge TE percentage data with centromere length
merged_data = pd.merge(total_percentage, centromere_data, on='Organism')

# Classify each organism by group based on the organism name
merged_data['Group'] = merged_data['Organism'].apply(
    lambda x: 'Cryptococcus' if 'Cryptococcus' in x else 'Kwoniella')

# Statistical Analysis
#--------------------------------------------------------------------------------
# Testing Normality of the data with Shapiro-Wilk test
_, p_value_te = stats.shapiro(merged_data['percentage'])
_, p_value_cen = stats.shapiro(merged_data['cen_length_bp'])

# Correlation Analysis (deciding on Pearson or Spearman based on normality test results)
if p_value_te < 0.05 or p_value_cen < 0.05:
    corr_coef, p_value = stats.spearmanr(merged_data['cen_length_bp'], merged_data['percentage'])
    method_used = "Spearman's rho"
else:
    corr_coef, p_value = stats.pearsonr(merged_data['cen_length_bp'], merged_data['percentage'])
    method_used = "Pearson's r"

# Visualization
#--------------------------------------------------------------------------------
plt.figure(figsize=(6, 6))
sns.scatterplot(x='cen_length_bp', y='percentage', hue='Group', data=merged_data,
                palette={'Cryptococcus': '#f7a38f', 'Kwoniella': '#9bd8f2'})
plt.xlabel('Centromere Length (bp)')
plt.ylabel('Percentage of TEs')
plt.title('Correlation between Centromere Length and TE Percentage')
plt.text(x=np.min(merged_data['cen_length_bp']), y=np.max(merged_data['percentage']),
         s=f'{method_used}: {corr_coef:.2f}, p-value: {p_value:.3f}',
         horizontalalignment='left', verticalalignment='top')

plt.rcParams['pdf.fonttype'] = 42
plt.savefig('TEs_vs_AvgCenLength_Correlation_with_Coefficient.pdf')
plt.show()

# Output Results and Save Data
#--------------------------------------------------------------------------------
print(f'Correlation Coefficient ({method_used}): {corr_coef}')
print(f'P-value: {p_value}')

# Prepare and save statistical results
stats_results = {
    'Correlation Coefficient': corr_coef,
    'P-value': p_value,
    'Normality Shapiro-Wilk test p-value (TE %)': p_value_te,
    'Normality Shapiro-Wilk test p-value (CEN length in bp)': p_value_cen
}
stats_df = pd.DataFrame([stats_results])
stats_df.to_csv('TEs_vs_AvgCenLength_Correlation_Stats.csv', index=False)

# Define a custom sorting order for organisms
organism_order = [
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

# Generate a sorting key based on the custom order. If an organism isn't found in the list, it gets a default index.
merged_data['Sort_Key'] = merged_data['Organism'].apply(lambda x: organism_order.index(x) if x in organism_order else len(organism_order))

# Sort the DataFrame by the custom order and drop the Sort_Key column as it's no longer needed after sorting
merged_data_sorted = merged_data.sort_values('Sort_Key').drop('Sort_Key', axis=1)

# Save the sorted merged data
merged_data_sorted.to_csv('TEs_vs_AvgCenLength_MergedData_Sorted.csv', index=False)

# Display file paths of saved files
stats_file_path = 'TEs_vs_AvgCenLength_Correlation_Stats.csv'
merged_data_file_path = 'TEs_vs_AvgCenLength_MergedData_Sorted.csv'
print(f'Saved Files:\n{stats_file_path}\n{merged_data_file_path}')