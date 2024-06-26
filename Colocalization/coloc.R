############################################################################
#
#        Colollization of each two traits among GWAS, pQTL and eQTL  
#          
#             coloc package (v5.1.0)
#
############################################################################

# The following code demonstrates the case for GWAS and pQTL

options(scipen = 1, digits = 2)
library(dplyr)
library(coloc)
library(gCMAP)
IndexSNP <- read.table(gzfile("Meta.GWAS.IndexSNP.txt.gz"), header=F, sep="\t"); colnames(IndexSNP) <- c("CHR","IndexSNP","IndexRSID")
GWAS <- read.table(gzfile("Meta.GWAS.txt.gz"), header=F, sep="\t");   colnames(GWAS) <- c("IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N")
pQTL <- read.table(gzfile("Kidney.pQTL.txt.gz"), header=F, sep="\t"); colnames(pQTL) <- c("RSID","Gene","pQTL_BETA","pQTL_SE","pQTL_PVAL","pQTL_N")
Sig.pQTL <- read.table(gzfile("Kidney.pQTL.Sig.txt.gz"), header=T, sep="\t")

### Significant pQTLs based on Protein level FDR(qvalue) < 0.05
pQTL_forPair <- pQTL
pQTL_forPair$RSID_ProteinID <- paste(pQTL_forPair$RSID, pQTL_forPair$Protein, sep="_")
pQTL_forPair <- subset(pQTL_forPair, RSID_ProteinID %in% Sig.pQTL$RSID_ProteinID)[,1:6]

####################################################################################################
##### Colollization between eGFR GWAS and kidney pQTL 
GWAS.pQTL.coloc <- data.frame(matrix(nrow = 0, ncol = 25))  ### for final results
colnames(GWAS.pQTL.coloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","Protein","pQTL_BETA","pQTL_SE","pQTL_PVAL","pQTL_N","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")

### Run each indexSNP
for (Index in 1:nrow(IndexSNP)){
	### Combine all data based on POS because all in the same chromsome
	pQTL_forCurrent <- pQTL_forPair
	
	CurrentIndexSNP <- as.character(IndexSNP[Index,2])
    CurrentIndex <- IndexSNP[Index,]
    CurrentIndex <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));

    CurrentIndex_GWAS_pQTL <- inner_join(CurrentIndex, pQTL_forPair, by = c("RSID" = "RSID"));
    
    if(nrow(CurrentIndex_GWAS_pQTL) == 0) {
    	CurrentIndex_GWAS_pQTL <- inner_join(CurrentIndex, pQTL, by = c("RSID" = "RSID"));
    	Min_PVAL <- min(CurrentIndex_GWAS_pQTL$pQTL_PVAL)
    	Min_PVAL_Protein <- subset(CurrentIndex_GWAS_pQTL, pQTL_PVAL == Min_PVAL)
    	Min_PVAL_Protein <- as.character(Min_PVAL_Protein[1,]$Protein)
    	pQTL_forCurrent <- subset(pQTL, Protein == Min_PVAL_Protein)
    }
    	
	### Combine
    CurrentIndex <- inner_join(CurrentIndex, pQTL_forCurrent, by = c("RSID" = "RSID"));
    CurrentIndex <- na.omit(CurrentIndex);
    
    Pair <- unique(CurrentIndex[,c("IndexRSID","Protein")]); nrow(Pair)   ### Sig pQTL only
     
    Number <- nrow(Pair)
    currentIndex.GWAS.pQTL.coloc <- data.frame(matrix(nrow = 0, ncol = 25))
    colnames(currentIndex.GWAS.pQTL.coloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","Protein","pQTL_BETA","pQTL_SE","pQTL_PVAL","pQTL_N","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
	
	### Run each pair
    for (i in 1:Number){
		Current.pQTL <- subset(pQTL, Protein == Pair[i,2])
		
		CurrentIndex <- IndexSNP[Index,]
    	Current.pair <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));
   		Current.pair <- inner_join(Current.pair, Current.pQTL, by = c("RSID" = "RSID"));
    	Current.pair <- na.omit(Current.pair);
    	Current.pair$GWAS_Var <- (Current.pair$GWAS_SE)^2
    	Current.pair$pQTL_Var <- (Current.pair$pQTL_SE)^2
    			
		### Only do coloc if more than 50 SNP available 
  		if (nrow(Current.pair) >= 50) {
  			### Conbine SNP_Protein as ID
  			Current.pair$CHR <- paste(Current.pair$IndexRSID, Current.pair$Protein, sep="|"); colnames(Current.pair)[1] = "Pair"
  			
  			### IndexingSNP
			Current.pair.IndexSNP <- subset(Current.pair, as.character(IndexRSID) == as.character(RSID))
			if (nrow(Current.pair.IndexSNP) >=1){
				Current.GWAS.pQTL.coloc <- Current.pair.IndexSNP[,1:19]
			}else{ ## Indexing SNP does not have at least one kind of three data (maybe pQTL), use the SNP with lowest p value instead
				Current.GWAS.pQTL.coloc <- Current.pair[which.min(Current.pair$GWAS_PVAL),1:19]
			}
			Current.GWAS.pQTL.coloc[,20:25] <- NA
			colnames(Current.GWAS.pQTL.coloc) <- colnames(currentIndex.GWAS.pQTL.coloc)
						
			### coloc calculation
  			Current.gwas <- Current.pair[,c("RSID","POS","GWAS_BETA","GWAS_Var","GWAS_PVAL","GWAS_N","MAF")]; colnames(Current.gwas) <- c("snp","position","beta","varbeta","pvalues","N","MAF"); Current.gwas$snp <- as.character(Current.gwas$snp)
  			Current.pQTL <- Current.pair[,c("RSID","POS","pQTL_BETA","pQTL_Var","pQTL_PVAL","pQTL_N","MAF")]; colnames(Current.pQTL) <- c("snp","position","beta","varbeta","pvalues","N","MAF"); Current.pQTL$snp <- as.character(Current.pQTL$snp)
						
			Current.gwas <- as.list(as.data.frame(Current.gwas)); Current.gwas[["type"]] <- "quant"; str(Current.gwas)
			Current.pQTL <- as.list(as.data.frame(Current.pQTL)); Current.pQTL[["type"]] <- "quant"; str(Current.pQTL)
  			coloc <- coloc.abf(dataset1 = Current.gwas, dataset2 = Current.pQTL)     ### used for further analysis (default parameters)
  			
  			### coloc.Result
  			Current.GWAS.pQTL.coloc[1,20:25] <- coloc$summary
			currentIndex.GWAS.pQTL.coloc <- rbind(currentIndex.GWAS.pQTL.coloc,Current.GWAS.pQTL.coloc)
		} ## if
	} ### Run each pair
	GWAS.pQTL.coloc <- rbind(GWAS.pQTL.coloc,currentIndex.GWAS.pQTL.coloc)  ## Add current indexSNP result to finial result
} ### Run each indexSNP

##### Output final result
write.table(GWAS.pQTL.coloc, "GWAS.pQTL.coloc.txt", quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")
