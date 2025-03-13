#!/usr/bin/env Rscript
################################################################################
# Colocalization analysis (GWAS and QTL) using coloc package (coloc.susie)
#
# This script performs colocalization analysis using coloc.susie,
# leveraging SuSiE regression for fine-mapping.
#
# Required input:
#   - GWAS summary statistics
#   - QTL summary statistics
#   - Reference genotype data for LD calculation (based on the 1000 Genomes Project European population)
#   - Pair information for loci and phenotypes
#
# Output will be stored in the specified results directory.
################################################################################

suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(bigreadr))
suppressMessages(library(arrow))
suppressMessages(library(coloc))
suppressMessages(library(qvalue))
suppressMessages(library(susieR))
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))

options(stringsAsFactors = FALSE)

# Define command-line options
args_list <- list(
  make_option("--outdir", type = "character", default = "./", help = "Output directory"),
  make_option("--prefix", type = "character", default = "coloc", help = "Output prefix"),
  make_option("--input", type = "character", default = NULL, help = "QTL statistics (QTLtools output)"),
  make_option("--ref", type = "character", default = NULL, help = "Reference genotype data for LD"),
  make_option("--gwas", type = "character", default = NULL, help = "GWAS summary statistics"),
  make_option("--pair", type = "character", default = NULL, help = "Pair information for loci and phenotypes"),
  make_option("--chr", type = "numeric", default = 1, help = "Chromosome of interest"),
  make_option("--N", type = "numeric", default = NULL, help = "Sample size for QTL analysis"),
  make_option("--bed", type = "character", default = NULL, help = "Phenotype matrix")
)

opt_parser <- OptionParser(option_list = args_list)
opt <- parse_args(opt_parser)

# Set parameters
outdir <- opt$outdir
QTL_n <- opt$N
input <- opt$input
prefix <- opt$prefix
i <- opt$chr

# Load GWAS data
gwas <- fread2(opt$gwas) %>%
  rename_with(~ c('chr', 'start', 'end', 'MarkerName', 'REF', 'ALT', 'Freq_ALT', 'MAF', 'BETA', 'SE', 'P', 'N'))
gwas$P <- as.numeric(gwas$P)

gwas <- gwas %>% mutate(variant_id = paste(chr, start, REF, ALT, sep = ":"))

# Load QTL data
stat <- fread2(input) %>%
  rename_with(~ sub("Beta_ALT", "slope", .)) %>%
  rename_with(~ sub("SE", "slope_se", .)) %>%
  rename_with(~ sub("AF_ALT", "af", .)) %>%
  mutate(variant_id = paste(chr_QTL, start_QTL, REF_QTL, ALT_QTL, sep = ":"))

# Load pair information
pairs <- readRDS(opt$pair)

# Get list of loci
loci <- unique(pairs$rsid)

for (locus in loci) {
  phenos <- unique(pairs %>% filter(rsid == locus) %>% pull(phenotype_id))
  
  for (pheno in phenos) {
    # Filter GWAS-QTL overlap
    tmp_coloc <- stat %>%
      filter(phenotype_id == pheno, !is.na(slope), !is.na(slope_se)) %>%
      inner_join(gwas, by = "variant_id")
    
    if (nrow(tmp_coloc) < 25) next
    
    # Calculate sdY from phenotype matrix
    pheno_matrix <- fread2(opt$bed) %>% rename(Phenotype_ID = phenotype_id)
    v_sdY <- sd(na.omit(pheno_matrix %>% filter(Phenotype_ID == pheno) %>% select(-c(1:6)) %>% unlist()))
    if (is.na(v_sdY)) next
    
    # Prepare QTL dataset
    tmp_QTL_df <- tmp_coloc %>%
      transmute(snp = variant_id, phenotype_id, af, beta = slope, varbeta = slope_se^2, MAF = pmin(af, 1 - af), N = QTL_n, sdY = v_sdY)
    tmp_QTL_df$type <- "quant"
    
    # Prepare GWAS dataset
    tmp_gwas_df <- tmp_coloc %>%
      transmute(snp = variant_id, af = Freq_ALT, beta = BETA, varbeta = SE^2, pvalues = P, N, MAF)
    tmp_gwas_df$type <- "quant"
    
    # Compute LD matrix
    tmp_output <- file.path(outdir, paste0("tmp_", prefix, "_", i, "_", locus, "_", pheno))
    write.table(tmp_coloc$variant_id, file = tmp_output, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    cmd <- paste('plink2 --bfile', paste0(opt$ref, 'chr', i), '--r2 square', '--extract', tmp_output, '--out', tmp_output)
    system(cmd)
    
    LD <- fread2(paste0(tmp_output, ".unphased.vcor2"))
    rownames(LD) <- colnames(LD) <- fread2(paste0(tmp_output, ".unphased.vcor2.vars"))
    LD <- as.matrix(LD)
    
    common_snps <- intersect(tmp_QTL_df$snp, tmp_gwas_df$snp)
    common_snps <- intersect(common_snps, colnames(LD))
    
    tmp_QTL_df <- filter(tmp_QTL_df, snp %in% common_snps)
    tmp_gwas_df <- filter(tmp_gwas_df, snp %in% common_snps)
    LD <- LD[common_snps, common_snps, drop = FALSE]
    
    if (length(common_snps) < 25) next
    
    # Run susie_rss
    S1 <- susie_rss(tmp_gwas_df$beta / sqrt(tmp_gwas_df$varbeta), LD, max(tmp_gwas_df$N), max_iter = 200)
    S2 <- susie_rss(tmp_QTL_df$beta / sqrt(tmp_QTL_df$varbeta), LD, max(tmp_QTL_df$N), max_iter = 200)
    
    susie.res <- coloc.susie(S1, S2)
    
    if (!is.null(susie.res$summary)) {
      output_file <- paste0(outdir, "/", prefix, "_coloc_susie_res_chr", i, "_", locus, "_", pheno, ".txt")
      write.table(susie.res$summary, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
  }
}
