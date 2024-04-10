
#!/bin/bash

# Script Metadata
#--------------------------------------------------------------------------------------------------------------
# Script name: CalculateGCContent_1kbWindows.sh
# Title: Calculate GC Content in 1kb Non-Overlapping Windows
# Author: Marco A. Coelho @ Heitman lab
# Date: 2023-04-06
# Description: This script calculates GC content in 1kb non-overlapping windows
#              for genomic sequences (*.fa files) using seqkit v2.0.0. The results are
#              consolidated into a single file with a header, ready for further analysis.
#
# Usage: Run this script in a directory containing .fa files. Ensure seqkit is installed
#        and the seqkit_env_v2.0.0 Conda environment is correctly set up.
#        Adjust the output_dir variable as necessary to specify output location.
#
# Requirements: Conda environment named seqkit_env_v2.0.0 with seqkit installel
#--------------------------------------------------------------------------------------------------------------

# Activate the required Conda environment with seqkit installed
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate seqkit_env_v2.0.0 || { echo "Failed to activate seqkit_env_v2.0.0"; exit 1; }

# Directory where output files will be stored (adjust as necessary)
output_dir="./GC_content_analysis"
mkdir -p "$output_dir"

# Loop through all fasta files (*.fa) in the current directory
for f in *.fa; do
    STRAIN=$(basename "$f" .scaffolds.fa)

    cat "$f" | \
    seqkit sliding -s 1000 -W 1000 | \
    seqkit fx2tab -n -g -i | \
    awk -v OFS='\t' -v strain="$STRAIN" '{print strain, $2}' > "$output_dir/${STRAIN}.GCcontent_1kb_window.txt"

    # Check if seqkit and awk commands succeeded
    if [ $? -ne 0 ]; then
        echo "Error processing file $f. Exiting."
        exit 1
    fi
done

# Add a header and combine all individual strain GC content files into a single file
{
    echo -e "Species\tGCperc"
    cat "$output_dir"/*.GCcontent_1kb_window.txt
} > "$output_dir/Cryptococcus_Kwoniella_GCcontent_1kb_window.txt" || { echo "Error combining files. Exiting."; exit 1; }

# Completion message
echo "GC content calculation and file consolidation completed successfully."

# Deactivate Conda environment
conda deactivate


