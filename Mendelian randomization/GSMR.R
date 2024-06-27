###################################################################################
#
#        Bidirectional GSMR analysis of SBP/DBP and eGFR
#          
#             gsmr package (v1.1.0)
#
###################################################################################

# The following code demonstrates the case for SBP and eGFR

library(gsmr)
library(dplyr)

setwd("/path/to/working_directory")
load("SBP.eGFR.clumped.RData")
gsmr_data <- SBP.eGFR.clumped
head(gsmr_data)

# Estimate LD correlation matrix using R
snp_coeff_id = scan("/path/to/SBP.eGFR.clumped.xmat.gz", what="", nlines=1)
snp_coeff = read.table("/path/to/SBP.eGFR.clumped.xmat.gz", header=F, skip=2)
snp_id = Reduce(intersect, list(gsmr_data$SNP, snp_coeff_id))
gsmr_data = gsmr_data[match(snp_id, gsmr_data$SNP),]
snp_order = match(snp_id, snp_coeff_id)
snp_coeff_id = snp_coeff_id[snp_order]
snp_coeff = snp_coeff[, snp_order]

# Calculate the LD correlation matrix
ldrho = cor(snp_coeff)

# Check the size of the correlation matrix and double-check if the order of the SNPs in the LD correlation matrix is consistent with that in the GWAS summary data
colnames(ldrho) = rownames(ldrho) = snp_coeff_id
dim(ldrho)
ldrho[1:5,1:5]

# Standardization
snpfreq = gsmr_data$a1_freq
bzx = gsmr_data$bzx
bzx_se = gsmr_data$bzx_se
bzx_n = gsmr_data$bzx_n
std_zx = std_effect(snpfreq, bzx, bzx_se, bzx_n)
gsmr_data$std_bzx = std_zx$b
gsmr_data$std_bzx_se = std_zx$se
head(gsmr_data)

# GSMR analysis
bzx = gsmr_data$std_bzx
bzx_se = gsmr_data$std_bzx_se
bzx_pval = gsmr_data$bzx_pval
bzy = gsmr_data$bzy
bzy_se = gsmr_data$bzy_se
bzy_pval = gsmr_data$bzy_pval
n_ref = 503
gwas_thresh = 5e-8
single_snp_heidi_thresh = 0.01
multi_snps_heidi_thresh = 0.01
nsnps_thresh = 10
heidi_outlier_flag = T
ld_r2_thresh = 0.05
ld_fdr_thresh = 0.05
gsmr2_beta = 0
gsmr_results = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)

cat("The estimated effect of the exposure on outcome: ", gsmr_results$bxy)
cat("Standard error of bxy: ", gsmr_results$bxy_se)
cat("P-value for bxy: ", gsmr_results$bxy_pval)
cat("Indexes of the SNPs used in the GSMR analysis: ", gsmr_results$used_index[1:5], "...")
cat("Number of SNPs with missing estimates in the summary data: ", length(gsmr_results$na_snps))
cat("Number of non-significant SNPs: ", length(gsmr_results$weak_snps))
cat("Number of SNPs in high LD (LD rsq >", ld_r2_thresh, "): ", length(gsmr_results$linkage_snps))
cat("Number of pleiotropic outliers: ", length(gsmr_results$pleio_snps))

# Bi-directional GSMR analysis
bi_gsmr_results = bi_gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snps_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)
cat("Effect of risk factor on disease: ", bi_gsmr_results$forward_bxy)

# Save
save(gsmr_results, bi_gsmr_results, file = "SBP.eGFR.clumped.GSMR.Results.RData")
