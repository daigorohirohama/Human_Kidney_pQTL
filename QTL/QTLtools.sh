##################################################################################
#
#        pQTL/eQTL identification by QTLtools  
#
##################################################################################

### The following code demonstrates the case for pQTL

##################################################################################
#        #1: cis permutation
##################################################################################

module load QTLtools/1.3.1-11

QTLtools cis \
--vcf /path/to/vcf/vcf_file.vcf.gz \
--bed /path/to/bed/bed_file.bed.gz \
--cov /path/to/covariate_file/cov_file.txt \
--window 1000000 \
--permute 1000 \
--std-err \
--out /path/to/output_directory/cis_permutation_pQTL.txt


###################################################################################
#        #2: cis Nominal        
###################################################################################

module load QTLtools/1.3.1-11

QTLtools cis \
--vcf /path/to/vcf/vcf_file.vcf.gz \
--bed /path/to/bed/bed_file.bed.gz \
--cov /path/to/covariate_file/cov_file.txt \
--window 1000000 \
--nominal 1 \
--std-err \
--out /path/to/output_directory/cis_nominal_pQTL.txt


###################################################################################
#        #3: trans Nominal        
###################################################################################

module load QTLtools/1.3.1-11

QTLtools trans \
    --vcf /path/to/vcf/vcf_file.vcf.gz \
    --bed /path/to/bed/bed_file.bed.gz \
    --cov /path/to/covariate_file/cov_file.txt \
    --nominal \
    --threshold 1e-5 \
    --out /path/to/output_directory/trans_nominal_pQTL.txt
