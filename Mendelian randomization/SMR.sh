############################################################################
#
#        SMR and HEIDI test for eGFR GWAS and kidney pQTL/eQTL  
#
#             SMR: Summary-data-based Mendelian Randomization    
#            HEIDI: Heterogeneity in dependent instruments analysis
#             SMR package (v1.3.1)
#
############################################################################

# Only variants with a FDR(qvalue) less than 0.05 at the protein or gene level were used as instruments
# The following code demonstrates the case for pQTL
# The SNP~protein pairs include all SNPs that meet the criterion of protein level FDR(qvalue) < 0.05

for SNP_protein in $(cat SNP_protein_list.txt); do
  ./Software/SMR/smr_Linux --bfile ./1000GP_Phase3/EUR_phase3_MAF05 \
  --gwas-summary /path/to/eGFR.GWAS.ma \
  --beqtl-summary /path/to/Kidney.pQTL \
  --diff-freq-prop 0.1 \
  --peqtl-heidi 0.05 \
  --extract-target-snp-probe ${SNP_protein} \
  --out /path/to/output_directory/eGFR.GWAS.Kidney.pQTL.${SNP_protein}.smr
done
