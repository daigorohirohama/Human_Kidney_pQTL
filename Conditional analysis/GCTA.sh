################################################################################
# GCTA-COJO fine mapping analysis
#
# This script performs conditional and joint analysis (COJO) using GCTA.
# The input for --cojo-slct consists only of pQTL lead variants.
#
# Required input:
#   - Genotype reference data (PLINK format)
#   - Summary statistics file for pQTL (lead variants only)
#   - Conditional SNP file for COJO conditional analysis (from COJO-slct results)
#
# Output:
#   - Independent lead variants selected via --cojo-slct
#   - Conditional analysis results via --cojo-cond
################################################################################

#!/bin/bash

# Set parameters
GCTA="/path/to/gcta64"  
INPUT_DIR="data/input"
OUTPUT_DIR="data/output"
BFILE="data/genotype_data"
WINDOW_SIZE=250  # 250kb
COND_FILE="$INPUT_DIR/COJO_cond_input.txt"  # Extracted from COJO-slct (.jma.cojo)

############################################################################
# Step 1: COJO-slct (Lead SNP selection using pQTL lead variants only)
############################################################################

for CHR in {1..22}; do
    $GCTA --bfile $BFILE \
          --chr $CHR \
          --cojo-file $INPUT_DIR/COJO_input_CHR_${CHR}.txt \
          --cojo-slct \
          --cojo-wind $WINDOW_SIZE \
          --out $OUTPUT_DIR/COJO_output_chr${CHR} &
done
wait

############################################################################
# Step 2: COJO-cond (Conditional analysis)
############################################################################

for CHR in {1..22}; do
    $GCTA --bfile $BFILE \
          --chr $CHR \
          --cojo-file $INPUT_DIR/COJO_input_CHR_${CHR}.txt \
          --cojo-cond $COND_FILE \
          --out $OUTPUT_DIR/COJO_cond_output_chr${CHR} &
done
wait

