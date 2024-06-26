############################################################################
#
#        Multiple trait colollization of GWAS, pQTL and eQTL  
#          
#             coloc package (v0.1.0)
#
############################################################################

options(scipen = 1, digits = 2)
library(dplyr)
library(moloc)
library(gCMAP)

IndexSNP <- read.table(gzfile("Meta.GWAS.IndexSNP.txt.gz"), header=F, sep="\t"); colnames(IndexSNP) <- c("CHR","IndexSNP","IndexRSID")
GWAS <- read.table(gzfile("Meta.GWAS.txt.gz"), header=F, sep="\t");   colnames(GWAS) <- c("IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N")
pQTL <- read.table(gzfile("Kidney.pQTL.txt.gz"), header=F, sep="\t"); colnames(pQTL) <- c("RSID","Protein","pQTL_BETA","pQTL_SE","pQTL_PVAL","pQTL_N")
eQTL <- read.table(gzfile("Kidney.eQTL.txt.gz"), header=F, sep="\t"); colnames(eQTL) <- c("RSID","Gene","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N")
Sig.pQTL <- read.table(gzfile("Kidney.pQTL.Sig.txt.gz"), header=T, sep="\t")
Sig.eQTL <- read.table(gzfile("Kidney.eQTL.Sig.txt.gz"), header=T, sep="\t")

### Significant pQTLs based on protein level FDR (q value) < 0.05
pQTL_forPair <- pQTL
pQTL_forPair$RSID_Protein <- paste(pQTL_forPair$RSID, pQTL_forPair$Protein, sep="_")
pQTL_forPair <- subset(pQTL_forPair, RSID_Protein %in% Sig.pQTL$RSID_Protein)[,1:6]

### Significant eQTLs based on gene level FDR (q value) < 0.05
eQTL_forPair <- eQTL
eQTL_forPair$RSID_GeneID <- paste(eQTL_forPair$RSID, eQTL_forPair$Gene, sep="_")
eQTL_forPair <- subset(eQTL_forPair, RSID_GeneID %in% Sig.eQTL$RSID_GeneID)[,1:6]

### Save result
res.moloc <- data.frame(matrix(nrow = 0, ncol = 54))  ### for final results

### Run each indexSNP
for (Index in 1:nrow(IndexSNP)){
	### Combine all data based on POS because all in the same chromsome
	pQTL_forCurrent <- pQTL_forPair
	eQTL_forCurrent <- eQTL_forPair
	
	CurrentIndexSNP <- as.character(IndexSNP[Index,2])
    CurrentIndex <- IndexSNP[Index,]
    CurrentIndex <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));
    
    CurrentIndex_GWAS_pQTL <- inner_join(CurrentIndex, pQTL_forPair, by = c("RSID" = "RSID"));
    CurrentIndex_GWAS_eQTL <- inner_join(CurrentIndex, eQTL_forPair, by = c("RSID" = "RSID"));
    
    if(nrow(CurrentIndex_GWAS_eQTL) == 0) {
    	CurrentIndex_GWAS_eQTL <- inner_join(CurrentIndex, eQTL, by = c("RSID" = "RSID"));
    	Min_PVAL <- min(CurrentIndex_GWAS_eQTL$eQTL_PVAL)
    	Min_PVAL_Gene <- subset(CurrentIndex_GWAS_eQTL, eQTL_PVAL == Min_PVAL)
    	Min_PVAL_Gene <- as.character(Min_PVAL_Gene[1,]$Gene)
    	eQTL_forCurrent <- subset(eQTL, Gene == Min_PVAL_Gene)
    }
    
    if(nrow(CurrentIndex_GWAS_pQTL) == 0) {
    	CurrentIndex_GWAS_pQTL <- inner_join(CurrentIndex, pQTL, by = c("RSID" = "RSID"));
    	Min_PVAL <- min(CurrentIndex_GWAS_pQTL$pQTL_PVAL)
    	Min_PVAL_Protein <- subset(CurrentIndex_GWAS_pQTL, pQTL_PVAL == Min_PVAL)
    	Min_PVAL_Protein <- as.character(Min_PVAL_Protein[1,]$Protein)
    	pQTL_forCurrent <- subset(pQTL, Protein == Min_PVAL_Protein)
    }
	
	### Combine
    CurrentIndex <- inner_join(CurrentIndex, pQTL_forCurrent, by = c("RSID" = "RSID"));
    CurrentIndex <- inner_join(CurrentIndex, eQTL_forCurrent, by = c("RSID" = "RSID"));
    CurrentIndex <- na.omit(CurrentIndex);
    
    Pair <- unique(CurrentIndex[,c("IndexRSID","Protein","Gene")]); nrow(Pair)   ### Sig pQTL and eQTL only
    Number <- nrow(Pair)
    currentIndex.res.moloc <- data.frame(matrix(nrow = 0, ncol = 54))
    colnames(currentIndex.res.moloc) <- c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","Protein","pQTL_BETA","pQTL_SE","pQTL_PVAL","pQTL_N","Gene","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N","nsnps","PPA_a","PPA_a_b","PPA_a_c","PPA_a_bc","PPA_a_b_c","PPA_b","PPA_b_c","PPA_ac_b","PPA_c","PPA_ab_c","PPA_ab","PPA_ac","PPA_bc","PPA_abc","PPA_zero","Coloc_ppas_a","Coloc_ppas_b","Coloc_ppas_ab","Coloc_ppas_c","Coloc_ppas_ac","Coloc_ppas_bc","Coloc_ppas_abc","Best_snp_coloc_a","Best_snp_coloc_b","Best_snp_coloc_ab","Best_snp_coloc_c","Best_snp_coloc_ac","Best_snp_coloc_bc","Best_snp_coloc_abc")
	
	### Run each pair
    for (i in 1:Number){
		Current.pQTL <- subset(pQTL, Protein == Pair[i,2])
		Current.eQTL <- subset(eQTL, Gene == Pair[i,3])
		
		CurrentIndex <- IndexSNP[Index,]
    	Current.GWAS.pQTL.eQTL <- inner_join(CurrentIndex, GWAS, by = c("IndexRSID" = "IndexRSID"));
    	Current.GWAS.pQTL.eQTL <- inner_join(Current.GWAS.pQTL.eQTL, Current.pQTL, by = c("RSID" = "RSID"));
   		Current.GWAS.pQTL.eQTL <- inner_join(Current.GWAS.pQTL.eQTL, Current.eQTL, by = c("RSID" = "RSID"));
    	Current.GWAS.pQTL.eQTL <- na.omit(Current.GWAS.pQTL.eQTL);
		
		### Only do moloc if more than 10 SNP available 
  		if (nrow(Current.GWAS.pQTL.eQTL) >= 50) {
  			### Conbine SNP_Protein_Gene as ID
  			Current.GWAS.pQTL.eQTL$CHR <- paste(Current.GWAS.pQTL.eQTL$IndexRSID, Current.GWAS.pQTL.eQTL$Protein, Current.GWAS.pQTL.eQTL$Gene, sep="|"); colnames(Current.GWAS.pQTL.eQTL)[1] = "Pair"
  			
  			### IndexingSNP
			Current.GWAS.pQTL.eQTL.IndexSNP <- subset(Current.GWAS.pQTL.eQTL, as.character(IndexRSID) == as.character(RSID))
			if (nrow(Current.GWAS.pQTL.eQTL.IndexSNP) >=1){
				Current.res.moloc <- Current.GWAS.pQTL.eQTL.IndexSNP
			}else{ ## Indexing SNP does not have at least one kind of three data (maybe eQTL), use the SNP with lowest p value instead
				Current.res.moloc <- Current.GWAS.pQTL.eQTL[which.min(Current.GWAS.pQTL.eQTL$GWAS_PVAL),]
			}
			Current.res.moloc[,25:54] <- NA
			colnames(Current.res.moloc) <- colnames(currentIndex.res.moloc)
			
			### coloc calculation
  			gwas <- Current.GWAS.pQTL.eQTL[,c("RSID","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","MAF")]; gwas[,7] <- "quant"; colnames(gwas) <- c("SNP","BETA","SE","PVAL","N","MAF","type");
  			pQTL <- Current.GWAS.pQTL.eQTL[,c("RSID","pQTL_BETA","pQTL_SE","pQTL_PVAL","pQTL_N","MAF")]; pQTL[,7] <- "quant"; colnames(pQTL) <- c("SNP","BETA","SE","PVAL","N","MAF","type");
  			eqtl <- Current.GWAS.pQTL.eQTL[,c("RSID","eQTL_BETA","eQTL_SE","eQTL_PVAL","eQTL_N","MAF")]; eqtl[,7] <- "quant"; colnames(eqtl) <- c("SNP","BETA","SE","PVAL","N","MAF","type");

  			Current.GWAS.pQTL.eQTL.moloc  <- list(gwas,pQTL,eqtl)
  			moloc <- moloc_test(Current.GWAS.pQTL.eQTL.moloc)       ### used for further analysis (default parameters)
  			
  			# moloc
  			Current.res.moloc[1,25] <- moloc$nsnps
  			Current.res.moloc[1,26:40] <- moloc$priors_lkl_ppa[,4]
  			Current.res.moloc[1,41:47] <- moloc$best_snp[,1]
			Current.res.moloc[1,48:54] <- as.vector(moloc$best_snp[,2])
			currentIndex.res.moloc <- rbind(currentIndex.res.moloc,Current.res.moloc)
		} ## if
	} ### Run each pair
	res.moloc <- rbind(res.moloc,currentIndex.res.moloc)  ## Add current indexSNP result to finial result
} ### Run each indexSNP

### To match previous code with Zscore
### https://rdrr.io/bioc/gCMAP/src/R/CMAPResults-accessors.R
res.moloc$pQTL_Zscore <- zScores(res.moloc$pQTL_PVAL, direction=res.moloc$pQTL_BETA); 
res.moloc$eQTL_Zscore <- zScores(res.moloc$eQTL_PVAL, direction=res.moloc$eQTL_BETA); 
res.moloc <- res.moloc[,c("Pair","IndexSNP","IndexRSID","RSID","SNP","POS","REF","ALT","F","MAF","GWAS_BETA","GWAS_SE","GWAS_PVAL","GWAS_N","Protein","pQTL_Zscore","pQTL_PVAL","pQTL_N","Gene","eQTL_Zscore","eQTL_PVAL","eQTL_N","pQTL_BETA","pQTL_SE","eQTL_BETA","eQTL_SE","nsnps","PPA_a","PPA_a_b","PPA_a_c","PPA_a_bc","PPA_a_b_c","PPA_b","PPA_b_c","PPA_ac_b","PPA_c","PPA_ab_c","PPA_ab","PPA_ac","PPA_bc","PPA_abc","PPA_zero","Coloc_ppas_a","Coloc_ppas_b","Coloc_ppas_ab","Coloc_ppas_c","Coloc_ppas_ac","Coloc_ppas_bc","Coloc_ppas_abc","Best_snp_coloc_a","Best_snp_coloc_b","Best_snp_coloc_ab","Best_snp_coloc_c","Best_snp_coloc_ac","Best_snp_coloc_bc","Best_snp_coloc_abc")]

##### Output final result
write.table(res.moloc, "res.moloc.txt", quote=FALSE, row.names=FALSE, col.names=TRUE,sep="\t")
