##Plots and result visualization
library(boa)
library(coda)

chain <- read.table("../out/Selection/bayescan/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_80anc_t.sel", colClasses ="numeric")
chain <- mcmc(chain,thin=10)
plot(chain)
summary(chain)
autocorr.diag(chain)
autocorr.plot(chain)
effectiveSize(chain)
geweke.diag(chain, frac1=0.1, frac2=0.5)

parameter="Fst2"
for (parameter in colnames(sel)[-1]) {
  plot(trace(sel[[parameter]]), xlab=parameter, main= paste(parameter, "posterior distribution"))
}

boa.hpd(sel[[parameter]], 0.05)


data <- read.table("../out/Selection/bayescan/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_t_fst.txt", colClasses = "numeric")
sel <- read.table("../out/Selection/bayescan/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_t.sel", colClasses ="numeric")
loci <- read.delim("../data/bed/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_thinn.bim", header = F)

source("BayeScan2.1/BayeScan2.1/R functions/plot_R.r")

qalpha <- 0.01
lentil.bayes  <- plot_bayescan("../out/Selection/bayescan/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_t_fst.txt", FDR= qalpha)
## alpha > 0 means positive or divergent selection, alpha < 0 indicates negative selection.
## I'm looking for divergent selection
outliers.bayes <- droplevels(loci[data$qval < qalpha,])
outliers.bayes[,4] <- as.numeric(as.character(outliers.bayes$V4))

######################################################
###Identify genes 
genes <- read.delim("../data/Lculinaris/Lens_culinaris_2.0.genes_description.txt")
genes$Start <- as.numeric(as.character(genes$Start))
genes$End <- as.numeric(as.character(genes$End))
genes$Chr <- as.character(genes$Chr)

genes.outl <- numeric(0)

for (i in 1:nrow(outliers.bayes)) {
  if (nrow(genes[genes$Chr == as.character(outliers.bayes[i,1]) & genes$Start <= outliers.bayes[i,4] & genes$End >= outliers.bayes[i,4],])>0) {
    a<- genes[genes$Chr == outliers.bayes[i,1] & genes$Start <= outliers.bayes[i,4] & genes$End >= outliers.bayes[i,4],]
    aa <- cbind(a, outliers.bayes[i,c(4)])
    genes.outl <- rbind(genes.outl, aa)
  }
}
genes.outl <- genes.outl[,-c(8,9)]

#Just genes
genes.id <- droplevels(unique(genes.outl[genes.outl$Type=="gene",-c(8:10)]))

write.table(genes.outl, file = "../out/Selection/bayescan/genes_outliers_details_bayes.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(genes.id, file = "../out/Selection/bayescan/genes_outliers_bayes.txt", quote = F, row.names = F, col.names = T, sep = "\t")










