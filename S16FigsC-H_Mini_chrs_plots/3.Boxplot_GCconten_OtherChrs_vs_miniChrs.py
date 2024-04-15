#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: GC Content comparison between mini-chrs. and other chromosomes in Kwoniella species
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-10-10
Description: Calculates GC content from GenBank sequences, classifies sequences by chromosome type into groups
             based on external lists of 'Mini-chromosomes' compares groups using the Mann-Whitney U test, and
             visualizes differences. Outputs include detailed statistical tests and GC content plots.
Requirements: BioPython, NumPy, SciPy, Matplotlib, Seaborn

Input files:
    - GenBank File: Genomic data in GenBank format.
    - Mini Chromosomes File: Mini-chromosome identifiers.

Output files:
    - GC Content Plots: <filename>.GCperc_comparison_boxplots.png and .pdf
    - GC Content Data: <filename>_gc_contents.tsv
    - Statistical Tests Results: <filename>_GCcontent_stat_tests.tsv
---------------------------------------------------------------------------------------------------------------------
"""

import sys
import os
from Bio import SeqIO
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

def read_chromosomes(file_path):
    """Returns a list of chromosome names read from a specified file."""
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]

def classify_contigs(gc_contents, mini_chromosome_ids):
    """Classifies contigs based on a list of mini-chromosome identifiers."""
    mini_gc_contents = [gc_contents[ctg] for ctg in mini_chromosome_ids if ctg in gc_contents]
    others_gc_contents = [gc for ctg, gc in gc_contents.items() if ctg not in mini_chromosome_ids]
    return mini_gc_contents, others_gc_contents

def calculate_gc_content(sequence):
    """Calculates and returns the GC content percentage of a given sequence."""
    return (float(sequence.count('G') + sequence.count('C')) * 100.0) / len(sequence)

def parse_genbank_file(genbank_file):
    """Parses a GenBank file and returns a dictionary of contig identifiers and their GC contents."""
    gc_content_dict = {}
    with open(genbank_file, "r") as file:
        for record in SeqIO.parse(file, "genbank"):
            gc_content = calculate_gc_content(record.seq)
            gc_content_dict[record.id] = gc_content
    return gc_content_dict

def perform_mann_whitney_u_test(group1, group2):
    """Performs a Mann-Whitney U test between two groups and returns the U statistic and p-value."""
    u_statistic, p_value = stats.mannwhitneyu(group1, group2)
    return {'U_statistic': u_statistic, 'p_value': p_value}

def save_statistical_results(test_results, filename):
    """Saves statistical test results to a TSV file."""
    with open(f"{filename}_GCcontent_stat_tests.tsv", "w") as file:
        for key, value in test_results.items():
            file.write(f"{key}\t{value}\n")

def plot_gc_content_comparison(group1, group2, p_value, filename):
    """Generates and saves a box plot comparing GC content between two groups."""
    data = [group2, group1]  # Keep this order to match the labels
    labels = ["Other chrs.", "Mini-chrs."]
    colors = ["#c5717e", "#83afd3"]  # Custom colors for each group

    fig, ax = plt.subplots(figsize=(4, 4))
    plt.rcParams['pdf.fonttype'] = 42

    # Box plot with custom colors
    sns.boxplot(data=data, ax=ax, palette=colors, width=0.5)
    
    # Individual data points (stripplot) with neutral color and transparency
    sns.stripplot(data=data, color="white", ax=ax, jitter=0.1, marker='o', edgecolor='black', linewidth=1.5, size=8, alpha=0.7)
   
    # Add a horizontal line connecting the two groups and display the p-value above the line
    y_max = max([max(lst) for lst in data]) + 0.5  # Get the y-coordinate slightly above the highest data point
    ax.plot([0, 1], [y_max, y_max], color="black", lw=1.0)
    ax.text(0.5, y_max + 0.01, f'p = {p_value:.4f}', ha='center', fontsize=9, va='bottom')

    ax.set_xticklabels(labels)
    ax.set_ylabel("GC Content (%)", fontsize=12)
    ax.set_title(filename.replace("_", " "), ha='center', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(f'{filename}.GCperc_comparison_boxplots.png', dpi=300)
    plt.savefig(f'{filename}.GCperc_comparison_boxplots.pdf', dpi=300)
    plt.show()

def save_gc_content_to_tsv(gc_contents, base_name):
    """Saves GC content data to a TSV file."""
    output_file = f"{base_name}_gc_contents.tsv"
    with open(output_file, "w") as file:
        file.write("Contig\tGC_Content\n")
        for ctg, gc in gc_contents.items():
            file.write(f"{ctg}\t{gc}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script_name.py <path_to_genbank_file> <path_to_mini_chrs_file>")
        sys.exit(1)

    genbank_file, mini_chrs_file = sys.argv[1:3]

    mini_chromosome_ids = read_chromosomes(mini_chrs_file)

    gc_contents = parse_genbank_file(genbank_file)
    mini_gc_contents, others_gc_contents = classify_contigs(gc_contents, mini_chromosome_ids)
    test_results = perform_mann_whitney_u_test(mini_gc_contents, others_gc_contents)
    save_statistical_results(test_results, os.path.splitext(genbank_file)[0])
    plot_gc_content_comparison(mini_gc_contents, others_gc_contents, test_results['p_value'], os.path.splitext(genbank_file)[0])
    save_gc_content_to_tsv(gc_contents, os.path.splitext(genbank_file)[0])

