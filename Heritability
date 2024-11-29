################################################################################
#   Heritability Analysis using GCTA-GREML
#
#   Description:
#   This script performs heritability analysis for 873 pQTLs identified in
#   kidney tissue. Variants with FDR < 0.05 were pre-selected for this analysis.
#   The heritability (V(G)/Vp) for each protein is calculated, and the average
#   heritability is reported.
#
#   Requirements:
#   - GCTA (version 1.94.1)
#   - PLINK
################################################################################

# Load necessary modules
module load gcta/1.94.1
module load plink

################################################################################
# Step 1: Perform linkage disequilibrium (LD) pruning
################################################################################
plink --bfile /path/to/genotype_file \
      --extract /path/to/filtered_variants_FDR_0.05.txt \
      --indep-pairwise 50 5 0.2 \
      --out /path/to/pruned_variants

################################################################################
# Step 2: Create genetic relationship matrix (GRM)
################################################################################
gcta64 --bfile /path/to/genotype_file \
       --extract /path/to/pruned_variants.prune.in \
       --make-grm \
       --out /path/to/kidney_grm

################################################################################
# Step 3: Perform heritability analysis and calculate the average V(G)/Vp
################################################################################
# Initialize: Prepare output file
output_file="/path/to/combined_heritability_results.txt"
echo -e "Protein_ID\tV(G)\tSE_V(G)\tV(e)\tSE_V(e)\tVp\tSE_Vp\tV(G)/Vp\tSE_V(G)/Vp\tPval" > $output_file

# Variables for averaging
sum_heritability=0
protein_count=873

# Loop through each protein and perform analysis
for i in $(seq 1 $protein_count); do
    result_file="/path/to/heritability_results_$i.hsq"
    
    # Run GCTA-GREML for each protein
    gcta64 --grm /path/to/kidney_grm \
           --pheno /path/to/phenotype_file.txt \
           --reml --mpheno $i \
           --out /path/to/heritability_results_$i
    
    # Extract values from the result file
    V_G=$(grep "V(G)" $result_file | awk '{print $2}')
    SE_V_G=$(grep "V(G)" $result_file | awk '{print $3}')
    V_e=$(grep "V(e)" $result_file | awk '{print $2}')
    SE_V_e=$(grep "V(e)" $result_file | awk '{print $3}')
    Vp=$(grep "Vp" $result_file | awk '{print $2}')
    SE_Vp=$(grep "Vp" $result_file | awk '{print $3}')
    VG_Vp=$(grep "V(G)/Vp" $result_file | awk '{print $2}')
    SE_VG_Vp=$(grep "V(G)/Vp" $result_file | awk '{print $3}')
    Pval=$(grep "Pval" $result_file | awk '{print $2}')
    
    # Add values to the output file
    echo -e "$i\t$V_G\t$SE_V_G\t$V_e\t$SE_V_e\t$Vp\t$SE_Vp\t$VG_Vp\t$SE_VG_Vp\t$Pval" >> $output_file

    # Accumulate heritability values for calculating the average
    sum_heritability=$(echo "$sum_heritability + $VG_Vp" | bc)
done

# Calculate and save the average heritability
average_heritability=$(echo "$sum_heritability / $protein_count" | bc -l)
average_output="/path/to/average_heritability_results.txt"
