library(ggplot2)
library(RIdeogram)
library(stringr)

#Load data
cnv <- read.delim("../data/xhmm/Lcul_CNV_chr.txt", sep = " ")
#Load annotation
gff3 <- read.delim("../data/Lculinaris/Lens_culinaris_2.0.genes.gff3", header = F)
gff3 <- droplevels(gff3[gff3$V1=="Lcu.2RBY.Chr1" | gff3$V1=="Lcu.2RBY.Chr2" | gff3$V1=="Lcu.2RBY.Chr3" | 
                            gff3$V1=="Lcu.2RBY.Chr4" | gff3$V1=="Lcu.2RBY.Chr5" | gff3$V1=="Lcu.2RBY.Chr6" | gff3$V1=="Lcu.2RBY.Chr7",])

bed <- read.delim("../data/Lculinaris/Lens_culinaris_2.0.genes.bed", header = F)

#################### GENES DENSITY ###########################
genes <- gff3[gff3$V3=="gene",]

genes$Rank <- cut(genes$V5, seq(0, 630000000, 1000000))
genes.dens <- as.data.frame(table(genes$V1, genes$Rank))
genes.dens <- genes.dens[genes.dens$Freq!=0,]
genes.dens <- cbind(genes.dens[,c(1,3)], str_split_fixed(genes.dens$Var2, pattern = ",", n=2))
colnames(genes.dens) <- c("Chr", "Value", "Start", "End")
genes.dens$Start <- as.numeric(as.character(gsub("\\(", "", paste(genes.dens$Start))))
genes.dens$End <- as.numeric(as.character(gsub("\\]", "", paste(genes.dens$End))))

#################### RESISTANCE GENES ###########################
resist <- read.delim("../data/ResistanceGenes/Lens_culinaris.RGA.info2.txt")
##Add start and end positions
resist.pos <- numeric(0)
for (i in 1:nrow(resist)) {
  a <- cbind(bed[bed$V4==as.character(resist[i,1]), c(1,2,3,4)])
  resist.pos <- rbind(resist.pos, a)
}
colnames(resist.pos) <- c("Chr", "Start", "End", "ID")
resist.pos <- droplevels(resist.pos[resist.pos$Chr=="Lcu.2RBY.Chr1" | resist.pos$Chr=="Lcu.2RBY.Chr2" | resist.pos$Chr=="Lcu.2RBY.Chr3" |
                         resist.pos$Chr=="Lcu.2RBY.Chr4" |resist.pos$Chr=="Lcu.2RBY.Chr5" |resist.pos$Chr=="Lcu.2RBY.Chr6" |resist.pos$Chr=="Lcu.2RBY.Chr7",])
resist.pos$Start <- as.numeric(as.character(resist.pos$Start))
resist.pos$End <- as.numeric(as.character(resist.pos$End))

##Resistance genes density
resist.pos$Rank <- cut(resist.pos$End, seq(0, 630000000, 1000000))
resist.dens <- as.data.frame(table(resist.pos$Chr, resist.pos$Rank))
resist.dens <- resist.dens[resist.dens$Freq!=0,]
resist.dens <- cbind(resist.dens[,c(1,3)], str_split_fixed(resist.dens$Var2, pattern = ",", n=2))
colnames(resist.dens) <- c("Chr", "Value", "Start", "End")
resist.dens$Start <- as.numeric(as.character(gsub("\\(", "", paste(resist.dens$Start))))
resist.dens$End <- as.numeric(as.character(gsub("\\]", "", paste(resist.dens$End))))

####################### CNVs ####################################
#CNVs density
variants <-  droplevels(unique(cnv[,c(2:8)]))
variants$Rank <- cut(variants$MID_BP, seq(0, 630000000, 1000000))
cnv.dens <- as.data.frame(table(variants$CHR, variants$Rank))
cnv.dens <- cnv.dens[cnv.dens$Freq!=0,]
cnv.dens <- cbind(cnv.dens[,c(1,3)], str_split_fixed(cnv.dens$Var2, pattern = ",", n=2))
colnames(cnv.dens) <- c("Chr", "Value", "Start", "End")
cnv.dens$Start <- as.numeric(as.character(gsub("\\(", "", paste(cnv.dens$Start))))
cnv.dens$End <- as.numeric(as.character(gsub("\\]", "", paste(cnv.dens$End))))

####################### DELETIONS ####################################
#Deletions density
del.dens <- as.data.frame(table(variants[variants$CNV=="DEL", 4], variants[variants$CNV=="DEL", 8]))
del.dens <- del.dens[del.dens$Freq!=0,]
del.dens <- cbind(del.dens[,c(1,3)], str_split_fixed(del.dens$Var2, pattern = ",", n=2))
colnames(del.dens) <- c("Chr", "Value", "Start", "End")
del.dens$Start <- as.numeric(as.character(gsub("\\(", "", paste(del.dens$Start))))
del.dens$End <- as.numeric(as.character(gsub("\\]", "", paste(del.dens$End))))

####################### DUPLICATIONS ####################################
#Deletions density
dup.dens <- as.data.frame(table(variants[variants$CNV=="DUP", 4], variants[variants$CNV=="DUP", 8]))
dup.dens <- dup.dens[dup.dens$Freq!=0,]
dup.dens <- cbind(dup.dens[,c(1,3)], str_split_fixed(dup.dens$Var2, pattern = ",", n=2))
colnames(dup.dens) <- c("Chr", "Value", "Start", "End")
dup.dens$Start <- as.numeric(as.character(gsub("\\(", "", paste(dup.dens$Start))))
dup.dens$End <- as.numeric(as.character(gsub("\\]", "", paste(dup.dens$End))))

###################### PLOT ####################################
## Karyotype
karyo <- numeric(0)
for (i in 1:7) {
  a <- cbind(levels(genes$V1)[i], 0, max(genes[genes$V1==levels(genes$V1)[i],5])+100000)
  karyo <- rbind(karyo, a)
}
colnames(karyo) <- c("Chr", "Start", "End")
karyo <- as.data.frame(karyo)
karyo$Start <- as.numeric(as.character(karyo$Start))
karyo$End <- as.numeric(as.character(karyo$End))
karyo$Chr <- as.character(karyo$Chr)

#### Plot ideogram
### CNVs 
ideogram(karyotype = karyo, overlaid = cnv.dens, label_type = "heatmap", 
         colorset1 = c("#b2d9b4", "#82cf87", "#006307"), Lx = 160, Ly = 190,
         output = "../out/xhmm/ideogram/ideo-cnv.svg")
convertSVG("../out/xhmm/ideogram/ideo-cnv.svg", file = "../out/xhmm/ideogram/ideo-cnv.png", device = "png")

### CNVs and resistance genes
ideogram(karyotype = karyo, overlaid = cnv.dens, label = resist.dens, label_type = "heatmap", 
         colorset1 = c("#b2d9b4", "#82cf87", "#006307"), colorset2 = c("#d6d6d6", "#867e91", "#454545"), Lx = 160, Ly = 190,
         output = "../out/xhmm/ideogram/ideo-cnv-resist.svg")
convertSVG("../out/xhmm/ideogram/ideo-cnv-resist.svg", file = "../out/xhmm/ideogram/ideo-cnv-resist.png", device = "png")

### CNVs and genes density
ideogram(karyotype = karyo, overlaid = cnv.dens, label = genes.dens, label_type = "heatmap", 
         colorset1 = c("#b2d9b4", "#82cf87", "#006307"), colorset2 = c("#e8e8e8", "#8c8c8c", "#121212"), Lx = 160, Ly = 190,
         output = "../out/xhmm/ideogram/ideo-cnv-genes.svg")
convertSVG("../out/xhmm/ideogram/ideo-cnv-genes.svg", file = "../out/xhmm/ideogram/ideo-cnv-genes.png", device = "png")

### Deletions
ideogram(karyotype = karyo, overlaid = del.dens, label_type = "heatmap", 
         colorset1 = c("#ffe9e8", "#cc4343", "#850000"), Lx = 160, Ly = 190,
         output = "../out/xhmm/ideogram/ideo-del.svg")
convertSVG("../out/xhmm/ideogram/ideo-del.svg", file = "../out/xhmm/ideogram/ideo-del.png", device = "png")

### Deletions and resistance genes
ideogram(karyotype = karyo, overlaid = del.dens, label = resist.dens, label_type = "heatmap", 
         colorset1 = c("#ffe9e8", "#cc4343", "#850000"), colorset2 = c("#d6d6d6", "#867e91", "#454545"), Lx = 160, Ly = 190,
         output = "../out/xhmm/ideogram/ideo-del-resist.svg")
convertSVG("../out/xhmm/ideogram/ideo-del-resist.svg", file = "../out/xhmm/ideogram/ideo-del-resist.png", device = "png")

### Deletions and genes density
ideogram(karyotype = karyo, overlaid = del.dens, label = genes.dens, label_type = "heatmap", 
         colorset1 = c("#ffe9e8", "#cc4343", "#850000"), colorset2 = c("#e8e8e8", "#8c8c8c", "#121212"), Lx = 160, Ly = 190,
         output = "../out/xhmm/ideogram/ideo-del-genes.svg")
convertSVG("../out/xhmm/ideogram/ideo-del-genes.svg", file = "../out/xhmm/ideogram/ideo-del-genes.png", device = "png")


### Duplication
ideogram(karyotype = karyo, overlaid = dup.dens, label_type = "heatmap", 
         colorset1 = c("#b8e3f2", "#43aacc", "#006485"), Lx = 160, Ly = 200,
         output = "../out/xhmm/ideogram/ideo-dup.svg")
convertSVG("../out/xhmm/ideogram/ideo-dup.svg", file = "../out/xhmm/ideogram/ideo-dup.png", device = "png")

### Duplication and resistance genes
ideogram(karyotype = karyo, overlaid = dup.dens, label = resist.dens, label_type = "heatmap", 
         colorset1 = c("#b8e3f2", "#43aacc", "#006485"), colorset2 = c("#d6d6d6", "#867e91", "#454545"), Lx = 160, Ly = 200,
         output = "../out/xhmm/ideogram/ideo-dup-resist.svg")
convertSVG("../out/xhmm/ideogram/ideo-dup-resist.svg", file = "../out/xhmm/ideogram/ideo-dup-resist.png", device = "png")

### Deletions and genes density
ideogram(karyotype = karyo, overlaid = dup.dens, label = genes.dens, label_type = "heatmap", 
         colorset1 = c("#b8e3f2", "#43aacc", "#006485"), colorset2 = c("#e8e8e8", "#8c8c8c", "#121212"), Lx = 160, Ly = 190,
         output = "../out/xhmm/ideogram/ideo-dup-genes.svg")
convertSVG("../out/xhmm/ideogram/ideo-dup-genes.svg", file = "../out/xhmm/ideogram/ideo-dup-genes.png", device = "png")

