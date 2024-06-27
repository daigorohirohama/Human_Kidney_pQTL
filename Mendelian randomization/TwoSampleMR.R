################################################################################
#
#        TwoSampleMR analysis for kidney pQTL/eQTL and plasma pQTL vs 38 GWAS data
#          
#             TwoSampleMR package (v0.5.6)
#
################################################################################

# The following code demonstrates the general case for TwoSampleMR analysis comparing kidney pQTLs against GWAS data

library(TwoSampleMR)

# Define input file paths
input_pqtl_gwas <- "path/to/input/pQTL_GWAS_data.txt"

# Read in the pQTL and GWAS data
data <- read.delim(input_pqtl_gwas, header=F, sep="\t")
colnames(data) <- c("SNP","EA_GWAS","BETA_GWAS","SE_GWAS","P_GWAS","SNP_pQTL","Protein","P_pQTL","BETA_pQTL","SE_pQTL","GeneSymbol","EA_pQTL")

# Align the effect alleles
data$BETA_pQTL_good <- ifelse(data$EA_pQTL == data$EA_GWAS, data$BETA_pQTL, -1 * data$BETA_pQTL)

# Initialize vectors to store results
BETA_MR <- NULL
SE_MR <- NULL
P_MR <- NULL

# Perform Mendelian Randomization using the Wald ratio method
for (i in 1:nrow(data)) {
  result <- mr_wald_ratio(data$BETA_pQTL_good[i], data$BETA_GWAS[i], data$SE_pQTL[i], data$SE_GWAS[i])
  BETA_MR <- c(BETA_MR, result$b)
  SE_MR <- c(SE_MR, result$se)
  P_MR <- c(P_MR, result$pval)
}

# Add MR results to the data
data$BETA_MR <- BETA_MR
data$SE_MR <- SE_MR
data$P_MR <- P_MR

# Define output file path
output_file <- "path/to/output/MR_results.txt"

# Save the results to a file
write.table(unique(data[order(data$P_MR), c(1:5,8,10,11,13:16)]), file=output_file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
