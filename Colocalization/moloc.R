#!/usr/bin/env Rscript
################################################################################
# Moloc analysis (GWAS, pQTL, eQTL)
#
# This script extracts candidate SNPs by overlapping GWAS, pQTL, and eQTL loci,
# and then performs moloc analysis.
#
# Required input:
#   - GWAS summary
#   - Significant GWAS loci
#   - pQTL summary
#   - Significant pQTL loci
#   - eQTL summary
#   - Significant eQTL loci
#
# Output will be stored in the results directory.
################################################################################

library(dplyr)
library(optparse)
library(bigreadr)
library(moloc)
library(future)

options(stringsAsFactors = FALSE)

# Define command-line options
args_list <- list(
  make_option("--gwas_all", type = "character", default = "data/gwas_summary.bed",
              help = "GWAS summary", metavar = "character"),
  make_option("--gwas_sig", type = "character", default = "data/gwas_sig.bed",
              help = "Significant GWAS loci", metavar = "character"),
  make_option("--qtl_dir", type = "character", default = "data/pqtl/",
              help = "Directory with pQTL summary files", metavar = "character"),
  make_option("--qtl_sig", type = "character", default = "data/pqtl_sig.bed",
              help = "Significant pQTL loci", metavar = "character"),
  make_option("--eqtl_dir", type = "character", default = "data/eqtl/",
              help = "Directory with eQTL summary files", metavar = "character"),
  make_option("--eqtl_sig", type = "character", default = "data/eqtl_sig.bed",
              help = "Significant eQTL loci", metavar = "character"),
  make_option("--outdir", type = "character", default = "results/",
              help = "Output directory", metavar = "character"),
  make_option("--prefix", type = "character", default = "moloc_res",
              help = "Prefix for output files", metavar = "character"),
  make_option("--threads", type = "numeric", default = 1,
              help = "Number of threads to use", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = args_list)
opt <- parse_args(opt_parser)

# Display options
cat("Using Options:\n")
cat("  --gwas_all:", opt$gwas_all, "\n")
cat("  --gwas_sig:", opt$gwas_sig, "\n")
cat("  --qtl_dir:", opt$qtl_dir, "\n")
cat("  --qtl_sig:", opt$qtl_sig, "\n")
cat("  --eqtl_dir:", opt$eqtl_dir, "\n")
cat("  --eqtl_sig:", opt$eqtl_sig, "\n")
cat("  --outdir:", opt$outdir, "\n")
cat("  --prefix:", opt$prefix, "\n")
cat("  --threads:", opt$threads, "\n\n")

# Check if files and directories exist
for (file in c(opt$gwas_all, opt$gwas_sig, opt$qtl_sig, opt$eqtl_sig)) {
  if (!file.exists(file)) stop(paste("ERROR: file not found:", file))
}
for (dir in c(opt$qtl_dir, opt$eqtl_dir)) {
  if (!file.exists(dir)) stop(paste("ERROR: directory not found:", dir))
}

if (!file.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# Setup parallel processing
future::plan("multicore", workers = opt$threads)
options(future.globals.maxSize = 1000 * 1024^2)

# Change working directory to output directory
setwd(opt$outdir)


# 1. Extract candidate SNPs ----------------------------------------------------
candidate_file <- paste0(opt$prefix, ".coloc.candidates.bed")
cmd <- paste("bedtools intersect -a", opt$qtl_sig, "-b", opt$gwas_sig, "-wa -wb >", candidate_file)
system(cmd)

# Read the candidate SNPs
candidates <- fread2(candidate_file) %>%
  rename(variant_id = V4) %>%
  distinct(variant_id, .keep_all = TRUE)

# Initialize result storage
res_moloc <- list()


# 2. Process each pQTL file ----------------------------------------------------
qtl_files <- list.files(opt$qtl_dir, pattern = "\\.txt$", full.names = TRUE)
qtl_files <- qtl_files[!grepl("SuSiE", qtl_files)]  # Exclude SuSiE-related files

for (qtl_file in qtl_files) {
  cat("\nProcessing pQTL file:", qtl_file, "\n")

  # Load pQTL data
  df_pQTL <- fread2(qtl_file) %>%
    filter(!is.na(slope)) %>%
    rename_with(~paste0("pQTL_", .), -variant_id) %>%
    distinct(variant_id, .keep_all = TRUE)

  # Load eQTL data
  eqtl_file <- paste0(opt$eqtl_dir, "/eQTL_", basename(qtl_file))
  if (!file.exists(eqtl_file)) next  # Skip if no matching eQTL file

  df_eQTL <- fread2(eqtl_file) %>%
    filter(!is.na(slope)) %>%
    rename_with(~paste0("eQTL_", .), -variant_id) %>%
    distinct(variant_id, .keep_all = TRUE)

  # Merge with GWAS summary
  tmp_gwas <- candidates %>% inner_join(fread2(opt$gwas_all), by = "variant_id")
  tmp_pQTL <- df_pQTL %>% filter(variant_id %in% tmp_gwas$variant_id)
  tmp_eQTL <- df_eQTL %>% filter(variant_id %in% tmp_gwas$variant_id)

  # Check if enough common SNPs exist
  common_snps <- Reduce(intersect, list(tmp_gwas$variant_id, tmp_pQTL$variant_id, tmp_eQTL$variant_id))
  if (length(common_snps) < 50) next  # Skip if fewer than 50 SNPs

  # Prepare data for moloc
  gwas <- tmp_gwas %>% transmute(SNP = variant_id, beta = BETA, varbeta = SE^2, pvalues = P, N, MAF, type = "quant")
  pqtl <- tmp_pQTL %>% transmute(SNP = variant_id, beta = pQTL_slope, varbeta = pQTL_slope_se^2, pvalues = pQTL_PVAL, N, MAF, type = "quant")
  eqtl <- tmp_eQTL %>% transmute(SNP = variant_id, beta = eQTL_slope, varbeta = eQTL_slope_se^2, pvalues = eQTL_PVAL, N, MAF, type = "quant")

  # Run moloc
  cat("\nRunning moloc for:", qtl_file, "\n")
  moloc_res <- moloc.abf(list(gwas = gwas, pQTL = pqtl, eQTL = eqtl))
  res_moloc[[qtl_file]] <- moloc_res$summary
}

# Save results
final_res <- do.call(rbind, res_moloc)
write.table(final_res, file = paste0(opt$prefix, ".moloc.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
