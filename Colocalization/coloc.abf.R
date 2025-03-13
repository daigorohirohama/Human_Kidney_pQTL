#!/usr/bin/env Rscript
################################################################################
# Colocalization analysis (GWAS and QTL) using coloc package (coloc.abf)
#
# This script extracts candidate SNPs by overlapping GWAS and QTL loci,
# and then performs colocalization analysis.
#
# Required input:
#   - GWAS summary
#   - Significant GWAS loci
#   - QTL summary
#   - Significant QTL loci
#
# Output will be stored in the results directory.
################################################################################

library(dplyr)
library(optparse)
library(bigreadr)
library(coloc)
library(qvalue)
library(future)

options(stringsAsFactors = FALSE)

# Define command-line options
args_list <- list(
  make_option("--gwas_all", type = "character", default = "data/gwas_summary.bed",
              help = "GWAS summary file (BED format)", metavar = "character"),
  make_option("--gwas_sig", type = "character", default = "data/gwas_significant_loci.bed",
              help = "GWAS significant loci file (BED format)", metavar = "character"),
  make_option("--qtl_dir", type = "character", default = "data/qtl/",
              help = "Directory containing QTL summary files", metavar = "character"),
  make_option("--qtl_sig", type = "character", default = "data/qtl_significant_loci.bed",
              help = "QTL significant loci file (BED format)", metavar = "character"),
  make_option("--outdir", type = "character", default = "results/",
              help = "Directory to store output files", metavar = "character"),
  make_option("--prefix", type = "character", default = "coloc_results",
              help = "Prefix for output files", metavar = "character"),
  make_option("--threads", type = "numeric", default = 1,
              help = "Number of threads for parallel processing", metavar = "numeric"),
  make_option("--mem", type = "numeric", default = 15,
              help = "Memory allocation (GB) for processing", metavar = "numeric")
)


opt_parser <- OptionParser(option_list = args_list)
opt <- parse_args(opt_parser)

# Display options
cat("Using Options:\n")
cat("  --gwas_all:", opt$gwas_all, "\n")
cat("  --gwas_sig:", opt$gwas_sig, "\n")
cat("  --qtl_dir:", opt$qtl_dir, "\n")
cat("  --qtl_sig:", opt$qtl_sig, "\n")
cat("  --outdir:", opt$outdir, "\n")
cat("  --prefix:", opt$prefix, "\n")
cat("  --threads:", opt$threads, "\n")
cat("  --mem:", opt$mem, "\n\n")

# Check if files and directories exist
if(!file.exists(opt$gwas_all)) stop(paste("ERROR: file not found:", opt$gwas_all))
if(!file.exists(opt$gwas_sig)) stop(paste("ERROR: file not found:", opt$gwas_sig))
if(!file.exists(opt$qtl_dir)) stop(paste("ERROR: directory not found:", opt$qtl_dir))
if(!file.exists(opt$qtl_sig)) stop(paste("ERROR: file not found:", opt$qtl_sig))

if(!file.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

# Setup parallel processing (if needed)
future::plan("multicore", workers = opt$threads)
options(future.globals.maxSize = 1000 * 1024^2)
options(future.seed = TRUE)

# Change working directory to output directory
setwd(opt$outdir)


# 1. Extract candidate SNPs ----------------------------------------------------
# Use bedtools to intersect QTL-significant regions and GWAS-significant loci.
# (Ensure bedtools is installed and in your PATH.)
candidate_file <- paste0(opt$prefix, ".coloc.candidates.bed")
cmd <- paste("bedtools intersect -a", opt$qtl_sig, "-b", opt$gwas_sig, "-wa -wb >", candidate_file)
system(cmd)

# Read the candidate file
candidates <- fread2(candidate_file)
colnames(candidates) <- c("chr", "pos", "pos2", "variant_id",
                          "CHR", "START", "END", "MarkerName", "locus")

# Initialize container for coloc results
res_coloc <- list()


# 2. Loop through QTL files ----------------------------------------------------
qtl_files <- list.files(opt$qtl_dir, pattern = "\\.txt$", full.names = TRUE)
# Exclude files that contain 'SuSiE' in the name (if any)
qtl_files <- qtl_files[!grepl("SuSiE", qtl_files)]

for(qtl_file in qtl_files) {
  
  # Create subdirectory for current QTL file results
  subdir <- tools::file_path_sans_ext(basename(qtl_file))
  subdir_path <- file.path(opt$outdir, subdir)
  if(!file.exists(subdir_path)) dir.create(subdir_path)
  setwd(subdir_path)
  
  # Load QTL summary stats; assume file has a column 'slope' and 'qval'
  qtl_df <- fread2(qtl_file) %>% filter(!is.na(slope))
  
  # Merge QTL data with candidate SNPs (based on variant_id)
  tmp_pairs <- qtl_df %>% inner_join(candidates, by = "variant_id") %>%
    filter(qval <= 0.05) %>% 
    transmute(CHR, START, END, locus, phenotype_id,
              pair = paste0(phenotype_id, "_", locus, "_", MarkerName)) %>%
    distinct(pair, .keep_all = TRUE)
  
  if(nrow(tmp_pairs) == 0) next
  
  # Process each candidate pair
  for(i in seq_len(nrow(tmp_pairs))) {
    tmp_bed <- tmp_pairs[i, ]
    tmp_prefix <- paste(tmp_bed$CHR, tmp_bed$pair, sep = "_")
    
    # Write temporary BED file defining the region
    bed_file <- paste0(tmp_prefix, ".bed")
    write.table(tmp_bed, file = bed_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    
    # Intersect GWAS summary stats with the region using bedtools
    gwas_tmp_file <- paste0(tmp_prefix, "_gwas.bed")
    cmd <- paste("bedtools intersect -a", opt$gwas_all, "-b", bed_file, "-wa -wb >", gwas_tmp_file)
    system(cmd)
    
    # Load GWAS summary stats from the intersected file
    tmp_gwas <- fread2(gwas_tmp_file)
    colnames(tmp_gwas) <- c("CHR", "START", "END", "MarkerName", "REF", "ALT",
                            "Freq_ALT", "MAF", "BETA", "SE", "P", "N",
                            "chr", "pos", "pos2", "locus", "phenotype_id", "pair")
    tmp_gwas <- tmp_gwas %>% 
      mutate(variant_id = paste(CHR, START, REF, ALT, sep = ":")) %>%
      distinct(variant_id, .keep_all = TRUE) %>%
      mutate(BETA = ifelse(Freq_ALT == MAF, BETA, -BETA),
             P = ifelse(P < 1e-320, 1e-320, P))
    
    # Retain only overlapping SNPs between GWAS and QTL datasets
    tmp_overlap <- tmp_gwas %>% inner_join(qtl_df, by = "variant_id")
    
    # Only proceed if at least 50 SNPs are available and some SNP is strongly significant
    if(nrow(tmp_overlap) < 50 ||
       nrow(tmp_overlap %>% filter(qval < 0.05 | P < 5e-08)) == 0) {
      next
    }
    
    # 3. Prepare data for coloc analysis ---------------------------------------
    # Build dataset list for GWAS (coloc.abf expects a list with required fields)
    gwas_list <- tmp_overlap %>% transmute(
      snp = variant_id,
      position = START,
      beta = BETA,
      varbeta = SE^2,
      pvalues = P,
      N,
      MAF
    ) %>% as.list
    gwas_list[["type"]] <- "quant"
    
    # Build dataset list for QTL
    qtl_list <- tmp_overlap %>% transmute(
      snp = variant_id,
      position = START,
      beta = slope,
      varbeta = slope_se^2,
      pvalues = pval_nominal,
      N = round(ma_count / af / 2),
      MAF = ifelse(af <= 0.5, af, 1 - af)
    ) %>% as.list
    qtl_list[["type"]] <- "quant"
    
    # 4. Run coloc ------------------------------------------------------------------
    coloc_res <- coloc.abf(dataset1 = gwas_list, dataset2 = qtl_list)
    
    # Sensitivity analysis: adjust p12 if the rule "H4 > 0.8" is met
    sens <- sensitivity(coloc_res, rule = "H4 > 0.8", doplot = FALSE)
    if(any(sens$pass == TRUE & sens$p12 >= 1e-05)) {
      new_p12 <- sens %>% filter(pass == TRUE, p12 >= 1e-05) %>% .$p12 %>% head(1)
      coloc_res <- coloc.abf(dataset1 = gwas_list, dataset2 = qtl_list, p12 = new_p12)
    }
    
    # Store the coloc summary along with metadata
    res_entry <- c(tmp_bed %>% transmute(chr = CHR, locus, phenotype_id, pair),
                   coloc_res$summary)
    res_coloc[[length(res_coloc) + 1]] <- as.data.frame(res_entry, stringsAsFactors = FALSE)
  } # End of candidate pair loop
  
  # Save intermediate results for current QTL file
  saveRDS(res_coloc, file = file.path(opt$outdir, paste0(opt$prefix, ".tmp.coloc.rds")))
  setwd(opt$outdir)
} # End of QTL file loop


# 5. Save final results ----------------------------------------------------------
final_res <- do.call(rbind, res_coloc)
saveRDS(final_res, file = paste0(opt$prefix, ".coloc.rds"))
write.table(final_res, file = paste0(opt$prefix, ".coloc.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
