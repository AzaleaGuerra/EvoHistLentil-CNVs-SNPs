#! /bin/bash

data="../data"
bayescan="BayeScan2.1/BayeScan2.1/binaries/BayeScan2.1_linux64bits"
pgdspider="PGDSpider_2.1.1.5/PGDSpider2-cli.jar"
ncores="15"
prefix="LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_80anc_thinn"

#Tranform from vcf to Bayescan format using PGDspider (file indicating population clustering is needed)
java -Xmx8g -jar ${pgdspider} -inputfile ${data}/vcf/${prefix}.recode.vcf -inputformat VCF -outputfile ${data}/Bayescan/${prefix} -outputformat GESTE_BAYE_SCAN -spid vcf2bayes.spid

##Bayescan 2.1
./${bayescan} -threads ${ncores} -snp ${data}/Bayescan/${prefix} -od ${data}/Bayescan/

#The default values in BayeScan for the MCMC algorithm parameters are quite conservative 
#and ensure good convergence in most cases. 
