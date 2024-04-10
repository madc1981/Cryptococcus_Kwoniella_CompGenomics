# Import necessary libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu, iqr, spearmanr

"""
---------------------------------------------------------------------------------------------------------------------
Script for performing statistical analysis and plotting of genomic features between groups (Cryptococcus and Kwoniella).
This includes comparison of groups using Mann-Whitney U tests, generating boxplots, and performing Spearman correlation analysis.
Author: Marco A. Coelho @ Heitman Lab
Date: 2023-05-22
Data input file: Cryptococcus_and_Kwoniella_genomic_features_agat.tsv
---------------------------------------------------------------------------------------------------------------------
"""

# Set the PDF font type for saving figures to ensure compatibility
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 9  # Default font size for plots

# Function to load data from a file
def load_data(filepath):
    return pd.read_csv(filepath, sep='\t')

# Load your data
data = load_data('Cryptococcus_and_Kwoniella_genomic_features_agat.tsv')

# Define colors for groups
colors = {"Cryptococcus": "#f9bfb2", "Kwoniella": "#bae4f6"}

# Function to perform the Mann-Whitney U test
def compare_groups_mann_whitney(data, group_col, num_col):
    group1, group2 = data[group_col].unique()
    data1 = data[data[group_col] == group1][num_col]
    data2 = data[data[group_col] == group2][num_col]
    stat, p = mannwhitneyu(data1, data2)
    return p

# Function to plot boxplots with specified font sizes
def plot_boxplot_with_points(data, group_col, num_col, colors, font_size=9):
    plt.figure(figsize=(6, 6))
    sns.boxplot(x=group_col, y=num_col, data=data, palette=colors, showfliers=True)
    sns.stripplot(x=group_col, y=num_col, data=data, color='black', edgecolor='black', size=6, jitter=True)
    plt.title(f'{num_col} by Group', fontsize=font_size)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)

# Perform the analysis for each numeric column
group_col = 'Group'
numeric_cols = ['Genome_size_mb', 'No_of_genes', 'Mean_intron_number_in_cds', 'Mean_intron_length_in_cds_in_bp']
results = []
summary_stats = []

for col in numeric_cols:
    p_value = compare_groups_mann_whitney(data, group_col, col)
    results.append((col, p_value))

    # Plot the boxplot with points and customize font size
    plot_boxplot_with_points(data, group_col, col, colors, font_size=9)

    # Calculate and annotate mean, median, and IQR for each group
    col_stats = []
    for group in data[group_col].unique():
        group_data = data[data[group_col] == group][col]
        mean = group_data.mean()
        median = group_data.median()
        iqr_value = iqr(group_data)
        col_stats.append([group, mean, median, iqr_value])

        # Add median and mean lines
        plt.axhline(mean, color='blue', linestyle='--', linewidth=1.0, 
                    xmin=(0.1 if group == data[group_col].unique()[0] else 0.6), 
                    xmax=(0.4 if group == data[group_col].unique()[0] else 0.9), 
                    label='Mean' if group == data[group_col].unique()[0] else "")
        
        # Annotate mean and median values
        plt.text(data[group_col].unique().tolist().index(group), mean, f'Mean: {mean:.2f}', color='blue', ha="center", va="bottom", fontsize=9)
        plt.text(data[group_col].unique().tolist().index(group), median, f'Median: {median:.2f}', color='green', ha="center", va="bottom", fontsize=9)

    summary_stats.append([col, col_stats])

    # Annotate Mann-Whitney p-value
    plt.annotate(f'Mann-Whitney p-value = {p_value:.4f}', xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9)
    
    # Save the figure as PDF
    plt.savefig(f'{col}_boxplot.pdf')
    plt.close()

# Save the results and summary statistics to files
results_df = pd.DataFrame(results, columns=['Column', 'P-Value'])
results_df.to_csv('Mann_whitney_test_results.csv', index=False)

# Save summary statistics to a CSV file
for stat in summary_stats:
    col_name, stats_data = stat
    stats_df = pd.DataFrame(stats_data, columns=['Group', 'Mean', 'Median', 'IQR'])
    stats_df.to_csv(f'{col_name}_stats.csv', index=False)

# Spearman correlation analysis and plotting
    
# Include dynamic axis limit calculations and plotting customization
def calculate_rounding_scale(value_range):
    #Dynamically determine the rounding scale based on the range of the data.
    if value_range < 10:
        return 0.5  # For very small ranges, round to nearest 0.5
    elif value_range < 50:
        return 1  # For small ranges, round to nearest 1
    elif value_range < 100:
        return 10  # For moderate ranges, round to nearest 10
    else:
        return 500  # For large ranges, round to nearest 100 or more dynamically

def calculate_axis_limits_with_margin(min_val, max_val, margin_factor=0.05):
    #Calculate axis limits with added margins and dynamically rounded.
    range_val = max_val - min_val
    margin = range_val * margin_factor
    round_to_nearest = calculate_rounding_scale(range_val + 2 * margin)  # Adjust rounding scale based on overall range with margin
    
    min_limit = np.floor((min_val - margin) / round_to_nearest) * round_to_nearest
    max_limit = np.ceil((max_val + margin) / round_to_nearest) * round_to_nearest
    
    return min_limit, max_limit

correlation_results = []

# Define blocks for Spearman correlation analysis, similar to the Mann-Whitney U test and boxplot sections
for col in numeric_cols[1:]:  # Excluding 'Genome_size_mb' for independent variable
    y_min_limit, y_max_limit = calculate_axis_limits_with_margin(data[col].min(), data[col].max(), margin_factor=0.05)

    for group in data[group_col].unique():
        group_data = data[data[group_col] == group]
        corr_coef, p_value = spearmanr(group_data['Genome_size_mb'], group_data[col])

        correlation_results.append((group, 'Genome_size_mb', col, corr_coef, p_value))

        # Plotting scatter plot with dynamically adjusted axis scales
        plt.figure(figsize=(6, 4))
        sns.scatterplot(x='Genome_size_mb', y=col, data=group_data, color=colors[group])
        x_min_limit, x_max_limit = calculate_axis_limits_with_margin(data['Genome_size_mb'].min(), data['Genome_size_mb'].max(), margin_factor=0.05)
        plt.xlim(x_min_limit, x_max_limit)
        plt.ylim(y_min_limit, y_max_limit)
        plt.title(f'Genome Size vs {col} in {group}', fontsize=9)
        plt.annotate(f'Spearman\'s ρ = {corr_coef:.2f}, p-value = {p_value:.4f}', xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9)
        plt.savefig(f'{group}_Genome_Size_vs_{col}_correlation.pdf')
        plt.close()

# Save Spearman's correlation results
correlation_df = pd.DataFrame(correlation_results, columns=['Group', 'Variable 1', 'Variable 2', 'Spearman\'s rho', 'P Value'])
correlation_df.to_csv('Spearman_correlation_analysis_results.csv', index=False)

print('Analysis complete. Figures and results saved.')
