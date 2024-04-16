#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: Orthogroup categorization and visualization
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-10-27
Description: This script categorizes orthogroups based on their presence across multiple species and visualizes the
             data using a stacked bar plot. It supports excluding specific species, exporting categorized data for
             selected species, and adjusting plot and file outputs based on provided prefixes.
Requirements: pandas, matplotlib, seaborn

Usage: python script_name.py input_file.tsv --exclude "species1,species2" --output_prefix "my_analysis" 
       --species_to_export "species3,species4"
Input files:
    - TSV file: Contains gene counts across various species.

Output files:
    - Gene count data: <output_prefix>_gene_counts.tsv
    - Stacked bar plot: <output_prefix>_stacked_bar_plot.pdf
    - Categorized orthogroups (optional): <output_prefix>_categorized_orthogroups_for_specific_species.tsv
---------------------------------------------------------------------------------------------------------------------
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns

# Categorize each orthogroup
def categorize_orthogroup(row):
    total_species = (row > 0).sum()
    total_genes = row.sum()
    
    if total_species == len(row):
        return "core"
    elif total_species == 1 and total_genes > 1:
        return "unique"
    else:
        return "shared"

def generate_stacked_bar_plot(args):
    # Load the input TSV file
    df = pd.read_csv(args.input_file, sep='\t', index_col=0, low_memory=False)

    # If there are species to exclude, drop them from the dataframe
    if args.exclude:
        exclude_species = args.exclude.split(',')
        df = df.drop(columns=exclude_species, errors='ignore')

    # Count number of genes for each species
    gene_counts_df = df.apply(lambda col: col.map(lambda x: len(str(x).split(',')) if pd.notna(x) else 0))

    # Categorize each orthogroup
    orthogroup_categories = gene_counts_df.apply(categorize_orthogroup, axis=1)
    aggregated_counts = gene_counts_df.groupby(orthogroup_categories).sum().transpose()
    
    # Grouping by categories and ordering of the stacked bar plot
    categories_order = ["core", "shared", "unique"]
    aggregated_counts = aggregated_counts[categories_order]
    
    # Reverse the order of the index before plotting to match the desire species ordering
    ordered_species_reversed = aggregated_counts.index.tolist()[::-1]

    # Plot settings
    plt.rcParams['pdf.fonttype'] = 42

    # Create a horizontal stacked bar plot with correct species ordering
    colormap = ["#024b7a", "#44a5c2", "#f15a29", "#ffae49" ]

    ax = aggregated_counts.reindex(ordered_species_reversed).plot(kind='barh', stacked=True, figsize=(8, 8), width=0.9, color=colormap)
    ax.xaxis.tick_top()                      # Move the x-axis to the top
    ax.xaxis.set_label_position('top')       # Move the x-axis label above the tick values
    ax.xaxis.labelpad = 12
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(True)
    ax.spines['top'].set_linewidth(1)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis="x", labelsize=10)
    ax.set_xlabel('Number of genes', fontsize=12)
    ax.legend(frameon=True, edgecolor='grey', fancybox=False)

    # Save gene count data to a TSV file
    output_file_path = f"{args.output_prefix}_gene_counts.tsv"
    aggregated_counts.to_csv(output_file_path, sep='\t')

    # Save plot to a PDF file
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_stacked_bar_plot.pdf",  bbox_inches='tight')
    plt.show()

    if args.species_to_export:
        species_to_export_list = args.species_to_export.split(',')
        export_categorized_orthogroups(args, species_to_export_list)

def export_categorized_orthogroups(args, species_list):
    df = pd.read_csv(args.input_file, sep='\t', index_col=0, low_memory=False)

    if args.exclude:
        exclude_species = args.exclude.split(',')
        df = df.drop(columns=exclude_species, errors='ignore')

    # Count genes for categorization
    gene_counts_df = df.apply(lambda col: col.map(lambda x: len(str(x).split(',')) if pd.notna(x) else 0))
    
    # Categorize each orthogroup
    orthogroup_categories = gene_counts_df.apply(categorize_orthogroup, axis=1)
    
    # Filter for the desired species
    filtered_df = df[species_list]
    
    # Filter out rows where all values for the specified species are zeros or NaN
    filtered_df = filtered_df[filtered_df[species_list].apply(lambda x: not all((y == 0 or pd.isna(y)) for y in x), axis=1)]

    # Add orthogroup category to the filtered dataframe
    filtered_df['Category'] = orthogroup_categories
    
    # Export the data
    output_file_path = f"{args.output_prefix}_categorized_orthogroups_for_specific_species.tsv"
    filtered_df.to_csv(output_file_path, sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate a stacked bar plot based on gene counts in core, shared, unique, and orphan categories.')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file (e.g. Orthogroups_renamed_sorted.tsv).')
    parser.add_argument('--exclude', type=str, help='Comma-separated list of species to exclude.', default=None)
    parser.add_argument('--output_prefix', type=str, help='Prefix for the output files.', default="output")
    parser.add_argument('--species_to_export', type=str, help='Comma-separated list of species to export categorized orthogroups/genes.', default=None)

    args = parser.parse_args()

    generate_stacked_bar_plot(args)
