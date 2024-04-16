#!/usr/bin/env python3
"""
Script Metadata
---------------------------------------------------------------------------------------------------------------------
Title: Script to process Orthogroup obtained from OrthoFinder for subsequent analysis
Author: Marco A. Coelho @ Heitman lab, Duke University
Date: 2023-10-27
Description: Processes orthogroup data by renaming columns based on a mapping list, reordering columns 
             based on a specified ordered list, and performing basic data integrity checks.
             The script ensures data is correctly formatted for downstream comparative genomic analysis.
Requirements: pandas

Input files:
            - Mapping_column_names.tsv: Contains mappings from original and final column header names (tab-separted)
            - Orthogroups.tsv: Contains orthogroup data to be processed (obtained from OrthoFinder, tab-separated)
            - Mapping_column_ordering.tsv: Contains mappings from original and final orders column orders (tab-separted).

Output files:
            - Orthogroups_renamed_sorted.tsv: Contains the processed orthogroup data with columns renamed and reordered.

Usage: python 0.rename_and_sort_OG_tables.py Mapping_column_names.tsv Orthogroups.tsv Mapping_column_ordering.tsv
---------------------------------------------------------------------------------------------------------------------
"""


import pandas as pd
import sys

def process_files(mapping_col_names, orthogroups_file, mapping_col_order):
    """Process and prepare orthogroup data for downstream analysis."""
    
    # Load the mapping from abbreviated to full column names
    column_mapping = pd.read_csv(mapping_col_names, sep="\t", header=None)
    column_mapping_dict = dict(zip(column_mapping[0], column_mapping[1]))

    # Read the orthogroup data and rename columns according to the mapping
    orthogroups = pd.read_csv(orthogroups_file, sep="\t")
    orthogroups_renamed = orthogroups.rename(columns=column_mapping_dict)

    # Load the final column order
    column_order = pd.read_csv(mapping_col_order, sep="\t", header=None)
    ordered_columns_names = column_order[1].tolist()
    final_ordered_columns = ['Orthogroup'] + [col for col in ordered_columns_names if col in orthogroups_renamed.columns]

    # Reorder the columns in the DataFrame
    orthogroups_reordered = orthogroups_renamed[final_ordered_columns]
    orthogroups_reordered.to_csv("Orthogroups_renamed_sorted.tsv", sep="\t", index=False)

    # Basic sanity checks
    assert orthogroups_reordered.columns[0] == 'Orthogroup'
    assert orthogroups_reordered.apply(lambda row: any(pd.notna(row[1:])), axis=1).all()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 0.rename_and_sort_OG_tables.py Mapping_column_names.tsv Orthogroups.tsv Mapping_column_ordering.tsv")
        sys.exit(1)
    
    mapping_col_names, orthogroups_file, mapping_col_order = sys.argv[1:4]
    process_files(mapping_col_names, orthogroups_file, mapping_col_order)