#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: Karyogram plotting of Gene genomic coverage for Kwoniella newhampshirensis CBS13917
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-10-11
Description: This script parses GenBank files to extract contig and gene lengths for a K. newhampshirensis where one
             chromosome is split into two contigs. These contigs are merged in the analysis. It visualizes the data
             in a karyogram, highlighting gene genomic coverage across contigs. Contigs are colored based on the
             ratio of gene length to contig length.
Requirements: BioPython, matplotlib, numpy

Input file:
    - GenBank file: Specifies the file containing genomic data in GenBank format.

Output files:
    - Karyogram plot: <filename>.gene_genomic_coverage.pdf, <filename>.gene_genomic_coverage.png
---------------------------------------------------------------------------------------------------------------------
"""

import sys
import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# Map ctgs to chr
contig_to_chromosome = {
    'ctg_1': 'chr_1', 'ctg_2': 'chr_2', 'ctg_3': 'chr_3', 'ctg_4': 'chr_4',
    'ctg_5': 'chr_5', 'ctg_6': 'chr_6', 'ctg_7': 'chr_7', 'ctg_8': 'chr_8',
    'ctg_9': 'chr_9', 'ctg_10': 'chr_11', 'ctg_11': 'chr_12', 'ctg_12': 'chr_13',
    'ctg_13': 'chr_14', 'ctg_14_15': 'chr_10', 'ctg_16': 'chr_15', 'ctg_17': 'chr_16'
}

def parse_genbank(file_name):
    """Parses a GenBank file, extracts contig and gene lengths, merges specific contigs, and renames contigs to chromosome identifiers.    
    This function handles specific genomic architecture by merging designated contigs (e.g., ctg_14 and ctg_15 into ctg_14_15)
    and then renames the contigs according to a predefined mapping to chromosome names to align with traditional chromosome nomenclature
    """
    contig_lengths = {}
    gene_lengths = {}
    with open(file_name, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contig_id = record.id
            contig_lengths[contig_id] = len(record)
            gene_lengths[contig_id] = sum(len(feature) for feature in record.features if feature.type == "gene")

    # Merging ctg_14 and ctg_15
    contig_lengths['ctg_14_15'] = contig_lengths['ctg_14'] + contig_lengths['ctg_15']
    gene_lengths['ctg_14_15'] = gene_lengths['ctg_14'] + gene_lengths['ctg_15']
    del contig_lengths['ctg_14'], contig_lengths['ctg_15']
    del gene_lengths['ctg_14'], gene_lengths['ctg_15']

    # Rename contigs to chromosomes
    contig_lengths = {contig_to_chromosome[k]: v for k, v in contig_lengths.items() if k in contig_to_chromosome}
    gene_lengths = {contig_to_chromosome[k]: v for k, v in gene_lengths.items() if k in contig_to_chromosome}

    return contig_lengths, gene_lengths


# Plotting the karyogram
def plot_karyogram(contig_lengths, gene_lengths):
    """Plots a karyogram showing contig lengths and gene genomic coverage"""
    # Extracts the base filename for use in plot titles and saving files
    filename = os.path.splitext(os.path.basename(file_name))[0]
    
	# Set the aesthetics for the plots
    plt.rcParams['pdf.fonttype'] = 42
    
    # Calculates gene genomic coverage as the ratio of total gene length to contig length for each contig
    gene_genomic_coverage = {contig: gene_lengths[contig] / contig_lengths[contig] for contig in contig_lengths}

    # Sorts contigs by length in descending order for consistent plotting
    sorted_contigs = sorted(contig_lengths.keys(), key=lambda x: contig_lengths[x], reverse=True)

    # Creates a list of coverage values for coloring the bars in the plot
    colors = [gene_genomic_coverage[contig] for contig in sorted_contigs]
    
    # Defines a color map and normalization for the colors based on coverage values
    cmap = plt.cm.RdBu_r
    norm = mcolors.Normalize(vmin=min(colors), vmax=max(colors))

    # Sets up the plot with a horizontal bar layout
    fig, ax = plt.subplots(figsize=(6, 5))
    y_pos = np.arange(len(sorted_contigs))

    # Maps normalized coverage values to colors using the colormap
    bar_colors = [cmap(norm(color)) for color in colors]

     # Plots horizontal bars with contig lengths converted to megabases
    ax.barh(y_pos, [contig_lengths[contig] / 1e6 for contig in sorted_contigs], color=bar_colors, zorder=2)

    # Sets y-axis labels to contig names and positions ticks
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sorted_contigs)

    # Labels the x-axis
    ax.set_xlabel('Contig Length (Mb)', fontsize=10)

    # Sets the title of the plot using the cleaned-up filename
    ax.set_title(filename.replace("_", " "), horizontalalignment='center', fontsize=12)

   # Adds a color bar to the plot to indicate the scale of genomic coverage
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical')
    cbar.set_label('Gene Genomic Coverage', size=10)  # Setting the font size for the label
    cbar.outline.set_visible(False)  # Hides the border of the color bar

    # Configures grid lines to appear behind plot elements
    ax.set_axisbelow(True)
    ax.grid(False)  # Disables grid lines
    ax.spines['top'].set_visible(False)  # Removes the top spine
    ax.spines['right'].set_visible(False)  # Removes the right spine
    ax.spines['left'].set_visible(False)  # Removes the left spine
    ax.spines['bottom'].set_visible(True)  # Removes the bottom spine

    # Inverts the y-axis so the largest contigs are at the top
    ax.invert_yaxis()

    # Adjusts layout to fit elements neatly
    plt.tight_layout()
    
	# Save the plot as PDF and PNG
    plt.savefig(f'{filename}.gene_genomic_coverage.pdf', bbox_inches='tight')
    plt.savefig(f'{filename}.gene_genomic_coverage.png', dpi=300, bbox_inches='tight')

    # Displays the plot
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide the GenBank file name as an argument.")
        sys.exit(1)
    file_name = sys.argv[1]
    contig_lengths, gene_lengths = parse_genbank(file_name)
    plot_karyogram(contig_lengths, gene_lengths)
