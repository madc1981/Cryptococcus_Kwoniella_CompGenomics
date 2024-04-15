#!/usr/bin/env python3

"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: Karyogram plotting of Gene genomic coverage
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-10-11
Description: This script parses GenBank files to extract contig and gene lengths and visualizes the data in a karyogram,
             highlighting gene genomic coverage across different contigs. The plot displays contig lengths and colors
             contigs based on the ratio of gene length to contig length.
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

def parse_genbank(file_name):
    """Parses a GenBank file to extract contig lengths and cumulative gene lengths for each contig"""
    contig_lengths = {}
    gene_lengths = {}
    with open(file_name, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contig_lengths[record.id] = len(record)
            gene_lengths[record.id] = sum(len(feature) for feature in record.features if feature.type == "gene")
    return contig_lengths, gene_lengths

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

    # Saves the plot to PDF and PNG with high resolution and tight bounding box
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
