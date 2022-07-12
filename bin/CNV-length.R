#library(xhmmScripts)
library(ggplot2)

#Load data. Only CNVs within chromosomes
cnv<- read.delim("../data/xhmm/Lcul_CNV_chr.txt", sep = "\t")
meta <- read.delim("../meta/LDP_meta.txt")

########### ACCORIDNG TO CNV ##########################
cnv$variant <- as.factor(paste(cnv$CNV, cnv$INTERVAL, sep = "-"))

## CNV Length
cnv$round.l <- round(cnv$KB, digits = -1)

#Keep sample with ancestryproportion > 0.80
cnv.80 <- droplevels(cnv[!is.na(cnv$GenGroup80),])

## Estimate median length according to genetics groups
median.gg <- numeric(0)
for (i in 1:8) {
  del.dup <- round(median(cnv[cnv$GenGroup==i,7]), 3)
  del.dup.90 <- round(median(cnv.90[cnv.90$GenGroup==i,7]), 3)
  del.dup.80 <- round(median(cnv.80[cnv.80$GenGroup==i,7]), 3)
  del <- round(median(cnv[cnv$GenGroup==i & cnv$CNV=="DEL",7]), 3)
  del.90 <- round(median(cnv.90[cnv.90$GenGroup==i & cnv.90$CNV=="DEL",7]), 3)
  del.80 <- round(median(cnv.80[cnv.80$GenGroup==i & cnv.80$CNV=="DEL",7]), 3)
  dup <- round(median(cnv[cnv$GenGroup==i & cnv$CNV=="DUP",7]), 3)
  dup.90 <- round(median(cnv.90[cnv.90$GenGroup==i & cnv.90$CNV=="DUP",7]), 3)
  dup.80 <- round(median(cnv.80[cnv.80$GenGroup==i & cnv.80$CNV=="DUP",7]), 3)
  m <- cbind(i, del.dup, del.dup.90, del.dup.80, del, del.90, del.80, dup, dup.90, dup.80)
  median.gg <- rbind(median.gg, m)
}

median.gg <- as.data.frame(rbind(median.gg, cbind("All", round(median(cnv$KB),3), 
                                                round(median(cnv.90$KB),3),
                                                round(median(cnv.80$KB),3),
                                                round(median(cnv[cnv$CNV=="DUP", 7]),3),
                                                round(median(cnv.90[cnv.90$CNV=="DUP", 7]),3),
                                                round(median(cnv.80[cnv.80$CNV=="DUP", 7]),3),
                                                round(median(cnv[cnv$CNV=="DEL", 7]),3),
                                                round(median(cnv.90[cnv.90$CNV=="DEL", 7]),3),
                                                round(median(cnv.80[cnv.80$CNV=="DEL", 7]),3))))
colnames(median.gg)[1] <- "GenGroup"
for (i in 2:ncol(median.gg)) {
  median.gg[,i] <- as.numeric(as.character(median.gg[,i]))
}

############ Perform a Kruskal-Wallis test ##########
kruskal.length.del <- kruskal.test(x=cnv[cnv$CNV=="DEL",7], 
                                   g= as.factor(cnv[cnv$CNV=="DEL",23]))
capture.output(pairwise.wilcox.test(x=cnv[cnv$CNV=="DEL",7], 
                                   g=as.factor(cnv[cnv$CNV=="DEL",23]), 
                                   p.adjust.method="holm"), file = "../out/xhmm/length/lenght_del_gg.posthoc")

kruskal.length.dup <- kruskal.test(x=cnv[cnv$CNV=="DUP",7], 
                                   g= as.factor(cnv[cnv$CNV=="DUP",23]))
capture.output(pairwise.wilcox.test(x=cnv[cnv$CNV=="DUP",7], 
                                   g=as.factor(cnv[cnv$CNV=="DUP",23]), 
                                   p.adjust.method="holm"), file = "../out/xhmm/length/lenght_dup_gg.posthoc")

## Violin plot/ Boxplot
ggplot(cnv, aes(x=CNV, y=KB, fill = CNV)) + 
  geom_boxplot(width=0.6, size = 0.5)+
  xlab(label = "CNV") +
  theme_bw(base_size = 25) 

ggplot(cnv, aes(x=as.factor(GenGroup), y=KB, fill = CNV)) + 
  #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), binwidth=5) +
  geom_boxplot(width=0.9, size = 0.5)+
  xlab(label = "Genetic groups")  + ylim(c(-250,3400)) +
  theme_bw(base_size = 25) 
