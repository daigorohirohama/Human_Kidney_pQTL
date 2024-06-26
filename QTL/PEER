#####################################################################################
#
#  PEER calculation for pQTL and eQTL
#     The following code demonstrates the case for pQTL
#     For eQTL, prepare the appropriate expression matrix and covariates, and calculate PEER factors using similar code
#
#####################################################################################

library(peer)
library(qtl)

########################################################################################
#  Prepare protein and covatiate file
########################################################################################
y <- read.csv("/path/to/protein_data_file.csv", header = TRUE, fileEncoding = "UTF-8-BOM", row.names = 1)  
 # The protein matrix has already been log2 transformed, followed by inverse normal transformation
y <- as.data.frame(t(y))  # Rows = ID, Columns = proteins
covs <- read.csv("/path/to/covariate_file.csv", header=TRUE, fileEncoding="UTF-8-BOM", stringsAsFactor=FALSE, row.names = 1)

#####################################################################################
#  Calculate PEER factors (Below shows the example case using K = 50)
#####################################################################################

Nmax_iterations <- 1000
K <- 50  

model <- PEER()
PEER_setNk(model, K)
PEER_setPhenoMean(model, as.matrix(y))
PEER_setPriorAlpha(model, 0.001, 0.1)
PEER_setPriorEps(model, 0.1, 10)
PEER_setNmax_iterations(model, Nmax_iterations)
PEER_setCovariates(model, as.matrix(covs))

PEER_update(model)

X <- PEER_getX(model)
W <- PEER_getW(model)
Alpha <- PEER_getAlpha(model)
Residuals <- PEER_getResiduals(model)
bounds <- PEER_getBounds(model)
vars <- PEER_getResidualVars(model)

colnames_X <- c("Age", "Gender", "Site", "eGFR", "PC1", "PC2", "PC3", "PC4", "PC5")
for (i in 1:K) {
  colnames_X <- append(colnames_X, paste("PEER", i, sep = ""))
}
colnames(X) <- colnames_X
rownames(X) <- rownames(covs)
rownames(Residuals) <- rownames(y)
colnames(Residuals) <- colnames(y)

output_prefix <- "PEER_results_K50"
output_dir <- "/path/to/output_directory/"

save(list = c("X", "W", "Alpha", "Residuals", "bounds", "vars"), 
     file = paste0(output_dir, output_prefix, ".RData"))
