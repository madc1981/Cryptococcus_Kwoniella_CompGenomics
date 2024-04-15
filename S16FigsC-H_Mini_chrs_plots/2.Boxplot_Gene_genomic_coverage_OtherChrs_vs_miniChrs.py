#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: Gene genomic coverage comparison between mini-chrs. and other chromosomes in Kwoniella species
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-10-12
Description: Calculates average gene length and gene Gene genomic coverage for each chromosome based on GFF3
             and chromosome size data. It classifies chromosomes into 'mini-chromosomes' and 'other chromosomes',
             performs Mann-Whitney U tests on these groups, and generates comparative boxplots.
Requirements: pandas, matplotlib, seaborn, scipy

Input files:
    - GFF3 file: Annotations for genes.
    - Chromosome sizes file: Contains sizes of chromosomes.
    - Mini Chromosomes File: List of mini-chromosome identifiers.

Output files:
    - Computed Metrics: <filename>_gene_metrics_output.tsv
    - Plots: <filename>.comparison_boxplots.png, <filename>.comparison_boxplots.pdf
---------------------------------------------------------------------------------------------------------------------
"""

import pandas as pd
import os
import sys
import re
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu


def natural_sort_key(s):
    """Helper function to sort strings containing numbers in a natural order."""
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

def read_chromosomes(file_path):
    """Reads a list of chromosome names from a file."""
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def compute_metrics(gff3_file, chr_sizes_file):
    """Computes average gene length and Gene genomic coverage for each chromosome"""
    
    # Load the GFF3 data
    gff3_data = pd.read_csv(gff3_file, sep='\t', comment='#', header=None,
                            names=['Contig', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Phase', 'Attributes'])
    
    # Load the chromosome sizes
    chr_sizes = pd.read_csv(chr_sizes_file, sep='\t', header=None, names=['Contig', 'Size'])

    # Filter the GFF3 data to retain only the gene entries
    gene_data = gff3_data[gff3_data['Type'] == 'gene'].copy()
    
    # Compute the length for each gene
    gene_data['Gene Length'] = gene_data['End'] - gene_data['Start'] + 1

    # Group by contig to compute the total gene length and number of genes for each contig
    grouped_gene_data = gene_data.groupby('Contig').agg(Total_Gene_Length=pd.NamedAgg(column='Gene Length', aggfunc='sum'),
                                                       Number_of_Genes=pd.NamedAgg(column='Type', aggfunc='count')).reset_index()
    
    # Merge the grouped data with chromosome sizes to get total size for each contig
    merged_data = pd.merge(grouped_gene_data, chr_sizes, on='Contig')
    
    # Compute the Average Gene Length for each contig (ratio of total gene length to total number of genes in contig)
    merged_data['Average Gene Length'] = merged_data['Total_Gene_Length'] / merged_data['Number_of_Genes']

    # Compute the Gene genomic coverage (ratio of total gene length to contig size) for each contig
    merged_data['Gene genomic coverage'] = merged_data['Total_Gene_Length'] / merged_data['Size']

    # Sort the data based on contig names in natural order
    merged_data = merged_data.sort_values(by='Contig', key=lambda x: x.map(natural_sort_key))
    return merged_data


def plot_and_statistical_tests(data, mini_chromosomes, output_filename):
    """Performs Mann-Whitney U tests and plots boxplots comparing two groups of chromosomes"""

    # Extract filename without extension
    filename = os.path.splitext(os.path.basename(gff3_file))[0]

    # Set the aesthetics for the plots
    sns.set(style="ticks")
    plt.figure(figsize=(8, 6))
    plt.rcParams['pdf.fonttype'] = 42
    plt.suptitle(f'{filename.replace("_", " ")}', fontsize=12, y=0.96)

    # Assign groups based on mini-chromosomes
    data['Group'] = data['Contig'].apply(lambda x: "Mini-chrs." if x in mini_chromosomes else "Other chrs.")
    group1 = data[data['Group'] == 'Mini-chrs.']
    group2 = data[data['Group'] == 'Other chrs.']

    # Collect test results
    test_results = []

    for idx, metric in enumerate(['Average Gene Length', 'Gene genomic coverage']):
        ax = plt.subplot(1, 2, idx + 1)
        sns.boxplot(x='Group', y=metric, data=data, palette=['#c5717e', '#83afd3'], width=0.5)
        sns.stripplot(x='Group', y=metric, data=data, color="white", jitter=0.1, marker='o', edgecolor='black', linewidth=1.5, size=8, alpha=0.7)
        u_stat, pvalue = mannwhitneyu(group1[metric], group2[metric])
        plt.title(f"P-value = {pvalue:.4f}", fontsize=11)
        plt.xlabel('')

        test_results.append({
            'Metric': metric,
            'U_Statistic': u_stat,
            'P_Value': pvalue
        })

        ax.tick_params(axis='x', which='major', length=5, direction='out')
        ax.tick_params(axis='y', which='major', length=5, direction='out')

    plt.tight_layout()
    plt.savefig(f'{filename}.GeneGenomicCov_comparison_boxplots.png', dpi=300)
    plt.savefig(f'{filename}.comparison_boxplots.pdf', dpi=300)
    plt.show()

    # Save test results to TSV
    pd.DataFrame(test_results).to_csv(output_filename, sep='\t', index=False)
    print(f"Statistical test results saved to {output_filename}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script_name.py <path_to_gff3_file> <path_to_chr_sizes_file> <path_to_mini_chrs_file>")
        sys.exit(1)

    gff3_file, chr_sizes_file, mini_chrs_file = sys.argv[1:4]
    mini_chromosomes = read_chromosomes(mini_chrs_file)
    results = compute_metrics(gff3_file, chr_sizes_file)
    results.to_csv(f"{os.path.splitext(gff3_file)[0]}_gene_metrics_output.tsv", sep="\t", index=False)
    print(f"Results saved to {os.path.splitext(gff3_file)[0]}_gene_metrics_output.tsv")

    # Specify filename for saving statistical results
    statistical_output_filename = f"{os.path.splitext(gff3_file)[0]}_GeneGenCov_stat_tests.tsv"
    plot_and_statistical_tests(results, mini_chromosomes, statistical_output_filename)
