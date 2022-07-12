library(stringr)
library(ggplot2)

#Load data
cnv <- read.delim("../data/CNVs/freq/cnv_freq.txt", sep = " ")
cnv.type <- read.delim("../data/CNVs/freq/cnv_freq_type.txt", sep = " ")

#Load annotation
genes <- read.delim("../data/Lculinaris/Lens_culinaris_2.0.genes_description.txt")
genes <- droplevels(genes[genes$Chr=="Lcu.2RBY.Chr1" | genes$Chr=="Lcu.2RBY.Chr2" | genes$Chr=="Lcu.2RBY.Chr3" | 
                          genes$Chr=="Lcu.2RBY.Chr4" | genes$Chr=="Lcu.2RBY.Chr5" | genes$Chr=="Lcu.2RBY.Chr6" | genes$Chr=="Lcu.2RBY.Chr7",])
genes$ID <- as.character(genes$ID)

bed <- read.delim("../data/Lculinaris/Lens_culinaris_2.0.genes.bed", header = F)
bed <- droplevels(bed[bed$V1=="Lcu.2RBY.Chr1" | bed$V1=="Lcu.2RBY.Chr2" | bed$V1=="Lcu.2RBY.Chr3" | 
                          bed$V1=="Lcu.2RBY.Chr4" | bed$V1=="Lcu.2RBY.Chr5" | bed$V1=="Lcu.2RBY.Chr6" | bed$V1=="Lcu.2RBY.Chr7",])

##Identify genes within CNVs
genes.cnv <- numeric(0)
for (i in 1:nrow(bed)) {
  for (j in 2:3) {
    if (nrow(cnv[cnv$CHR==bed[i,1] & cnv$Start<=bed[i,j] & cnv$End>=bed[i,j],])>0) {
      a<- cbind(cnv[cnv$CHR==bed[i,1] & cnv$Start<=bed[i,j] & cnv$End>=bed[i,j],],
                GeneID=bed[i,4])
      genes.cnv <- rbind(genes.cnv, a)
    }
  }
}
genes.cnv <- unique(genes.cnv)

#Bind gene description
genes.cnv.an <- numeric(0)
for (i in 1:nrow(genes.cnv)) {
  a <- cbind(genes.cnv[i,], genes[genes$ID==as.character(genes.cnv[i,7]), c(2, 6, 7)])
  genes.cnv.an <- rbind(genes.cnv.an, a)
}

#### List of genes' IDs (No frequencies)
gene.id <- droplevels(unique(genes.cnv[,7]))
write.table(gene.id, file = "../data/CNVs/genes/GenesIDs_CNV.txt", quote = F, col.names = F, row.names = F)

#### Count the presence of an affected gene by a CNV within a genetic group.
#A gene could be affected by different CNV with different frequencies at a given pop
counts.genes.gg <- numeric(0)
for (i in gene.id){
  for (j in 1:8) {
    a <- cbind(GenGroup=j, GeneID=i, Count=sum(genes.cnv.an[genes.cnv.an$GeneID==i & genes.cnv.an$GenGroup==j, 6]))
    counts.genes.gg <- rbind(counts.genes.gg, a)
  }
}
counts.genes.gg <- as.data.frame(counts.genes.gg)
counts.genes.gg[,3] <- as.numeric(as.character(counts.genes.gg[,3]))
counts.genes.gg <- counts.genes.gg[counts.genes.gg$Count > 0,]
write.table(counts.genes.gg, file = "../data/CNVs/genes/GenesCount_gg.txt", quote = F, col.names = T, row.names = F)

genes.id.10 <- droplevels(unique(counts.genes.gg[counts.genes.gg$Count>10, ]))
write.table(genes.id.10, file = "../data/CNVs/genes/GenesIDs_mac10.txt", quote = F, col.names = T, row.names = F)

##############################################################################
##Repete including the type of CNV (del or dup)
##Identify genes within CNVs
genes.cnv.type <- numeric(0)
for (i in 1:nrow(bed)) {
  for (j in 2:3) {
    if (nrow(cnv.type[cnv.type$CHR==bed[i,1] & cnv.type$Start<=bed[i,j] & cnv.type$End>=bed[i,j],])>0) {
      a<- cbind(cnv.type[cnv.type$CHR==bed[i,1] & cnv.type$Start<=bed[i,j] & cnv.type$End>=bed[i,j],],
                GeneID=bed[i,4])
      genes.cnv.type <- rbind(genes.cnv.type, a)
    }
  }
}
genes.cnv.type <- unique(genes.cnv.type)

#Bind gene description
genes.cnv.type.an <- numeric(0)
for (i in 1:nrow(genes.cnv.type)) {
  a <- cbind(genes.cnv.type[i,], genes[genes$ID==as.character(genes.cnv.type[i,8]), c(2, 6, 7)])
  genes.cnv.type.an <- rbind(genes.cnv.type.an, a)
}
write.table(genes.cnv.type.an, file = "../data/CNVs/genes/genes_gg_type.txt", quote = F, row.names = F, col.names = T, sep = "\t")
