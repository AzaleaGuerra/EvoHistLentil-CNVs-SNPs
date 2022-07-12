#!/usr/bin/env Rscript

###PCAdapt
library(pcadapt)
library(qvalue)
library(ggplot2)

#################################################
############ All accessions #####################
#Load data
lens <- read.pcadapt("../data/bed/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_thinn.bed", type = "bed")
meta <- read.delim("../meta/LDP_meta_SNPs.txt")

#PCA to evaluate the eigen vector that will be used
pca.lens <- pcadapt(lens, K= 20, min.maf = 0.0)

#plot explained variation by each eigen vector
jpeg("../out/Selection/pcadapt/var_explained.jpg", width = 900, height = 900)
plot(pca.lens, option= "screeplot")
dev.off()

jpeg("../out/Selection/pcadapt/pca.jpg", width = 900, height = 900)
plot(pca.lens, option = "scores", pop = meta$GenGroup)
dev.off()

#Repeat PCAdat considering only the first four eigenvector
pcadapt.lens <- pcadapt(lens, K= 7, min.maf = 0.0)
str(pcadapt.lens)

#Bonferroni correction
padj <- p.adjust(pcadapt.lens$pvalues,method="bonferroni")

#SNP position and SNP
pos.snps<- read.delim("../data/bed/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_thinn.bim", sep = "\t", header = FALSE)
pos.snps<- cbind(pos.snps[,c(1, 4)], paste(pos.snps[,1], pos.snps[,4], sep = "-"), 
                 pcadapt.lens$pvalues, padj)
colnames(pos.snps) <- c("CHR", "Pos", "SNP", "pval", "padj")

#Identify outliers
#Error rate 0.01/ #SNPs
alfa <- 0.05
lens.outliers <- na.omit(pos.snps[pos.snps$padj<alfa,])
nrow(lens.outliers)

write.table(lens.outliers, "../out/Selection/pcadapt/pcadapt_outliers", row.names = FALSE, quote = FALSE)

#Plots
jpeg("../out/Selection/pcadapt/p_dist-3.jpg", width = 900, height = 900)
hist(pcadapt.lens$pvalues, xlab = "p-values", main = NULL, breaks = 100, col = "orange")
dev.off()

jpeg("../out/Selection/pcadapt/manhattan-3.jpg", width = 900, height = 900)
plot(pcadapt.lens, option= "manhattan")+ 
  geom_point(data = lens.outliers, aes(x = as.numeric(rownames(lens.outliers)), y = -log10(pval)), color="red")
dev.off()

jpeg("../out/Selection/pcadapt/qqplot-3.jpg", width = 900, height = 900)
plot(pcadapt.lens, option= "qqplot")
  #geom_hline(yintercept = -log10(alfa), color="blue")
dev.off()


###################################################
############ Accessions > 0.80 ancestry ###########
#Load data
lens80 <- read.pcadapt("../data/bed/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_80anc_thinn.bed", type = "bed")
meta <- read.delim("../meta/LDP_meta_SNPs.txt")
meta80 <- meta[!is.na(meta$GenGroup80),]

#PCA to evaluate the eigen vector that will be used
pca.lens80 <- pcadapt(lens80, K= 20, min.maf = 0.0)

#plot explained variation by each eigen vector
jpeg("../out/Selection/pcadapt/var_explained80.jpg", width = 900, height = 900)
plot(pca.lens80, option= "screeplot")
dev.off()

jpeg("../out/Selection/pcadapt/pca80.jpg", width = 900, height = 900)
plot(pca.lens80, option = "scores", pop = meta80$GenGroup)
dev.off()

#Repeat PCAdat considering only the first three eigenvector
pcadapt.lens80 <- pcadapt(lens80, K= 6, min.maf = 0.0)
str(pcadapt.lens80)

#Error rate 0.01
alfa <- 0.05

#Bonferroni correction
padj80 <- p.adjust(pcadapt.lens80$pvalues,method="bonferroni")

#SNP position and SNP
pos.snps80<- read.delim("../data/bed/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV_80anc_thinn.bim", sep = "\t", header = FALSE)
pos.snps80<- cbind(pos.snps80[,c(1, 4)], paste(pos.snps80[,1], pos.snps80[,4], sep = "-"), 
                 pcadapt.lens80$pvalues, padj80)
colnames(pos.snps80) <- c("CHR", "Pos", "SNP", "pval", "padj")

#identify outliers
lens80.outliers <- na.omit(pos.snps80[pos.snps80$padj<alfa,])
nrow(lens80.outliers)

write.table(lens80.outliers, "../out/Selection/pcadapt/pcadapt_outliers80", row.names = FALSE, quote = FALSE)

#plots
jpeg("../out/Selection/pcadapt/p_dist80-3.jpg", width = 900, height = 900)
hist(pcadapt.lens80$pvalues, xlab = "p-values", main = NULL, breaks = 100, col = "orange")
dev.off()

jpeg("../out/Selection/pcadapt/manhattan80-3.jpg", width = 900, height = 900)
plot(pcadapt.lens80, option= "manhattan") + 
  geom_point(data = lens80.outliers, aes(x = as.numeric(rownames(lens80.outliers)), y = -log10(pval)), color="red")
dev.off()

jpeg("../out/Selection/pcadapt/qqplot80-3.jpg", width = 900, height = 900)
plot(pcadapt.lens80, option= "qqplot") 
  #geom_hline(yintercept = -log10(alfa), color="blue")
dev.off()


###Shared outlier snps 
lens.outliers <- read.delim("../out/Selection/pcadapt/pcadapt_outliers", sep = " ")
lens80.outliers <- read.delim("../out/Selection/pcadapt/pcadapt_outliers80", sep = " ")

library(VennDiagram)

venn.diagram(list(lens.outliers$SNP, lens80.outliers$SNP), 
             category.names = c("All accessions" , "Accessions > 0.80 ancestry"),
             filename = "../out/Selection/pcadapt/venn_pcadapt_SNPs.png",
             output = TRUE ,
             imagetype="png" ,
             resolution = 600,
             compression = "lzw",
             lwd = 1,
             lty = "blank",
             fill = c("#5BACED", "#F5AC2E"),
             cex = 1,
             fontfamily = "times",
             cat.cex = .7,
             #cat.default.pos = "outer",
             cat.pos = c(-15, 15),
             #cat.dist = c(0.055, 0.055, 0.055),
             cat.fontfamily = "times"
             )

