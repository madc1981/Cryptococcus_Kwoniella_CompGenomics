#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: Comparative Analysis of Transposable Element Content in Cryptococcus and Kwoniella Species
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-05-09
Description: This script analyzes the differences in transposable element (TE) content between Cryptococcus and 
Kwoniella species. It filters input data to exclude certain TE categories, calculates total TE percentages for each 
species, classifies species into Cryptococcus or Kwoniella groups, and performs statistical tests to compare TE 
content between these groups. The results are visualized in a boxplot, annotated with statistical test results, and 
saved along with a statistical summary.

Requirements: pandas, matplotlib, seaborn, scipy
Input file: Cryptococcus_Kwoniella_EarlGrey_TEs_merged.txt
Output files: TEs_Cryptococcus_vs_Kwoniella.pdf
              TEs_Cryptococcus_vs_Kwoniella_Stats.csv
              TEs_Cryptococcus_vs_Kwoniella_Table.csv
---------------------------------------------------------------------------------------------------------------------
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

# Load the data
file_path = 'Cryptococcus_Kwoniella_EarlGrey_TEs_merged.txt'
data = pd.read_csv(file_path, sep='\t')

# Filter out TE categories
boxplot_data = data[~data['tclassif'].isin(['Non-Repeat', 'Penelope', 'SINE'])]

# Calculating the total percentage of transposons for each organism
total_transposons = boxplot_data.groupby('Organism')['percentage'].sum().reset_index()

# Adding a column to classify each organism as Cryptococcus or Kwoniella
total_transposons['Group'] = total_transposons['Organism'].apply(lambda x: 'Cryptococcus' if 'Cryptococcus' in x else 'Kwoniella')

# Adjusting the 'percentage' values to reflect true percentages
total_transposons['percentage'] *= 100

# Organism order
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

# Create a sorting key based on the organism order
total_transposons['Sort_Key'] = total_transposons['Organism'].apply(lambda x: organism_order.index(x) if x in organism_order else len(organism_order))

# Sort the DataFrame by the custom order and drop the Sort_Key column
total_transposons_sorted = total_transposons.sort_values('Sort_Key').drop('Sort_Key', axis=1)

# Save the sorted raw data used for the boxplot
total_transposons_sorted.to_csv('TEs_Cryptococcus_vs_Kwoniella_Table.csv', index=False)

# Colors for each group
cryptococcus_color = "#f9bfb2"
kwoniella_color = "#bae4f6"

# Creating the boxplot with mean value and data points
plt.figure(figsize=(4, 8))
sns.boxplot(x='Group', y='percentage', data=total_transposons, 
            showmeans=True, meanline=True, 
            palette={'Cryptococcus': cryptococcus_color, 'Kwoniella': kwoniella_color})
sns.stripplot(x='Group', y='percentage', data=total_transposons, 
              color='grey', edgecolor='black', size=10, jitter=True)

# Statistical Tests
cryptococcus_transposons = total_transposons[total_transposons['Group'] == 'Cryptococcus']['percentage']
kwoniella_transposons = total_transposons[total_transposons['Group'] == 'Kwoniella']['percentage']

# Shapiro-Wilk Test for normality and Mann-Whitney U Test
p_cryptococcus = stats.shapiro(cryptococcus_transposons).pvalue
p_kwoniella = stats.shapiro(kwoniella_transposons).pvalue
u_stat, u_p_value = stats.mannwhitneyu(cryptococcus_transposons, kwoniella_transposons)

# Plot annotations with word-wrapped title
plt.title('Percentage of Transposons in\n Cryptococcus and Kwoniella', pad=10)
plt.xlabel('Group')
plt.ylabel('Percentage of TEs (%)')

# Adjusting subplot parameters and saving the plot as PDF
plt.subplots_adjust(top=0.92, bottom=0.06, left=0.157, right=0.969, hspace=0.2, wspace=0.2)
plt.rcParams['pdf.fonttype'] = 42
plt.savefig('TEs_Cryptococcus_vs_Kwoniella.pdf')
plt.show()

# Preparing statistical results for CSV
min_cryptococcus = cryptococcus_transposons.min()
max_cryptococcus = cryptococcus_transposons.max()
min_kwoniella = kwoniella_transposons.min()
max_kwoniella = kwoniella_transposons.max()

min_cryptococcus_species = total_transposons[(total_transposons['Group'] == 'Cryptococcus') & 
                                             (total_transposons['percentage'] == min_cryptococcus)]['Organism'].iloc[0]
max_cryptococcus_species = total_transposons[(total_transposons['Group'] == 'Cryptococcus') & 
                                             (total_transposons['percentage'] == max_cryptococcus)]['Organism'].iloc[0]
min_kwoniella_species = total_transposons[(total_transposons['Group'] == 'Kwoniella') & 
                                          (total_transposons['percentage'] == min_kwoniella)]['Organism'].iloc[0]
max_kwoniella_species = total_transposons[(total_transposons['Group'] == 'Kwoniella') & 
                                          (total_transposons['percentage'] == max_kwoniella)]['Organism'].iloc[0]

stats_results = {
    'Statistic': ['Mean Cryptococcus', 'Mean Kwoniella', 'Median Cryptococcus', 'Median Kwoniella',
                  'Min Cryptococcus', 'Min Kwoniella', 'Max Cryptococcus', 'Max Kwoniella',
                  'Shapiro-Wilk Cryptococcus', 'Shapiro-Wilk Kwoniella', 'Mann-Whitney U'],
    'Value': [cryptococcus_transposons.mean(), kwoniella_transposons.mean(), 
              cryptococcus_transposons.median(), kwoniella_transposons.median(),
              min_cryptococcus, min_kwoniella, max_cryptococcus, max_kwoniella,
              None, None, u_stat],
    'p-value': [None, None, None, None, None, None, None, None,
                p_cryptococcus, p_kwoniella, u_p_value],
    'Species': [None, None, None, None, min_cryptococcus_species, min_kwoniella_species, 
                max_cryptococcus_species, max_kwoniella_species, None, None, None]
}

stats_df = pd.DataFrame(stats_results)
stats_df.to_csv('TEs_Cryptococcus_vs_Kwoniella_Stats.csv', index=False)

print('Analysis complete. Files saved:')
print('TEs_Cryptococcus_vs_Kwoniella.pdf')
print('TEs_Cryptococcus_vs_Kwoniella_Stats.csv')
print('TEs_Cryptococcus_vs_Kwoniella_Table.csv')