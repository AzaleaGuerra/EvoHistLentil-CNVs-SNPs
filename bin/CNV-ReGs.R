library(stringr)
library(ggplot2)
library(RColorBrewer)

#Load data
cnv <- read.delim("../data/CNVs/freq/cnv_freq.txt", sep = " ")
cnv.type <- read.delim("../data/CNVs/freq/cnv_freq_type.txt", sep = " ")

#Load annotation
bed <- read.delim("../data/Lculinaris/Lens_culinaris_2.0.genes.bed", header = F)

genes <- read.delim("../data/Lculinaris/Lens_culinaris_2.0.genes_description.txt")
genes <- droplevels(genes[genes$Chr=="Lcu.2RBY.Chr1" | genes$Chr=="Lcu.2RBY.Chr2" | genes$Chr=="Lcu.2RBY.Chr3" | 
                            genes$Chr=="Lcu.2RBY.Chr4" | genes$Chr=="Lcu.2RBY.Chr5" | genes$Chr=="Lcu.2RBY.Chr6" | genes$Chr=="Lcu.2RBY.Chr7",])
genes$ID <- as.character(genes$ID)

#Resistance genes 
resist <- read.delim("../data/ResistanceGenes/Lens_culinaris.RGA.info2.txt")

##Add start and end positions
resist.pos <- numeric(0)
for (i in 1:nrow(resist)) {
  a <- cbind(bed[bed$V4==as.character(resist[i,1]), c(1,2,3,4)])
  resist.pos <- rbind(resist.pos, a)
}
colnames(resist.pos) <- c("Chr", "StartResGen", "EndResGen", "IDResGen")
resist.pos <- droplevels(resist.pos[resist.pos$Chr=="Lcu.2RBY.Chr1" | resist.pos$Chr=="Lcu.2RBY.Chr2" | resist.pos$Chr=="Lcu.2RBY.Chr3" |
                                      resist.pos$Chr=="Lcu.2RBY.Chr4" |resist.pos$Chr=="Lcu.2RBY.Chr5" |resist.pos$Chr=="Lcu.2RBY.Chr6" |resist.pos$Chr=="Lcu.2RBY.Chr7",])
resist.pos$StartResGen <- as.numeric(as.character(resist.pos$StartResGen))
resist.pos$EndResGen <- as.numeric(as.character(resist.pos$EndResGen))
resist.pos$KBResGen <- (resist.pos$End - resist.pos$StartResGen)/1000

## Find resistance genes within CNVs (genes are smaller than CNVs)
cnv.res <- numeric(0)
for (i in 1:nrow(resist.pos)) {
  for (j in 2:3) {
    if (nrow(cnv[cnv$CHR==resist.pos[i,1] & cnv$Start<=resist.pos[i,j] & cnv$End>=resist.pos[i,j],])>0) {
      a<- cbind(cnv[cnv$CHR==resist.pos[i,1] & cnv$Start<=resist.pos[i,j] & cnv$End>=resist.pos[i,j],],
                resist.pos[i,-c(1,5)])
      cnv.res <- rbind(cnv.res, a)
    }
  }
}

cnv.res <- unique(cnv.res)

#Assosciate funtional information
cnv.res.ann <- numeric(0)
for (i in 1:nrow(cnv.res)) {
  a <- cbind(cnv.res[i,], genes[genes$ID==as.character(cnv.res[i,9]), -c(10:14)])
  cnv.res.ann <- rbind(cnv.res.ann, a)
}
cnv.res.ann <- unique(cnv.res.ann)
write.table(cnv.res.ann, file = "../data/CNVs/RGA/CNV_ResistGenes.txt", quote = F, col.names = T, row.names = F, sep = "\t")

#### Count the presence of an affected gene by a CNV within a genetic group.
#A gene could be affected by different CNV with different frequencies at a given pop
counts.genes.gg <- numeric(0)
for (i in gene.id){
  for (j in 1:8) {
    a <- cbind(GenGroup=j, GeneID=i, Count=sum(cnv.res[cnv.res$IDResGen==i & cnv.res$GenGroup==j, 6]))
    counts.genes.gg <- rbind(counts.genes.gg, a)
  }
}
counts.genes.gg <- as.data.frame(counts.genes.gg)
counts.genes.gg[,3] <- as.numeric(as.character(counts.genes.gg[,3]))
counts.genes.gg <- counts.genes.gg[counts.genes.gg$Count > 0,]

#### List of resistance genes (No frequencies)
cnv.res <- droplevels(unique(cnv.res))
cnv.id <- droplevels(unique(cnv.res[,c(2)]))
del.id <- droplevels(unique(cnv.res[cnv.res$CNV=="DEL", c(1,2)]))
dup.id <- droplevels(unique(cnv.res[cnv.res$CNV=="DUP", c(1,2)]))
cnv.res.id <- droplevels(unique(cnv.res$IDResGen))
#write.table(cnv.res, file = "../out/xhmm/Resistance_genes/cnv-resistGenes.txt", quote = F, col.names = T, row.names = F)

