################################################################################
#
#        Hyprcoloc analysis for kidney pQTL/eQTL and plasma pQTL compared against 38 GWAS data
#          
#             Hyprcoloc implementation
#
################################################################################

# The following code demonstrates the general case for Hyprcoloc analysis comparing kidney pQTL against GWAS data

# Define input and output file paths
input_gwas="/path/to/input/gwas_data.tsv.bgz"
input_kidney_pqtl="/path/to/input/kidney_pqtl_data.txt.gz"
output_dir="/path/to/output/coloc/results"

# Define the trait labels
trait_labels="GWAS Kidney_pQTL"

# Run Hyprcoloc
/path/to/TargetLookup \
    --file ${input_gwas} ${input_kidney_pqtl} \
    --coloc \
    --labels ${trait_labels} \
    --build 37 37 \
    --dir ${output_dir}
