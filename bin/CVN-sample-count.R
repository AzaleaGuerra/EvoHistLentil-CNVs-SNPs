library("ggplot2")

#Load data. Only CNVs within chromosomes
cnv<- read.delim("../data/xhmm/Lcul_CNV_chr.txt", sep = "\t")
meta <- read.delim("../meta/LDP_meta.txt")

##############################################################
#Number of CNV according to sample
samp.c <- as.data.frame(table(cnv$SAMPLE, cnv$CNV))
#Bind metadata
samp.counts <- numeric(0)
for (i in 1:nrow(meta)) {
  a <- samp.c[samp.c$Var1 == as.character(meta[i,1]),]
  b <- cbind(a, meta[i,])
  samp.counts <- droplevels(rbind(samp.counts, b))
}

#Keep sample with ancestryproportion > 0.80
samp.counts.80 <- droplevels(samp.counts[!is.na(samp.counts$GenGroup80),])

#Get mean count according to genetic group
means.gg <- numeric(0)
for (i in 1:8) {
  del.dup <- round(mean(samp.counts[samp.counts$GenGroup==i,3]*2), 3)
  del.dup.80 <- round(mean(samp.counts.80[samp.counts.80$GenGroup==i,3]*2), 3)
  del <- round(mean(samp.counts[samp.counts$GenGroup==i & samp.counts$Var2=="DEL",3]), 3)
  del.80 <- round(mean(samp.counts.80[samp.counts.80$GenGroup==i & samp.counts.80$Var2=="DEL",3]), 3)
  dup <- round(mean(samp.counts[samp.counts$GenGroup==i & samp.counts$Var2=="DUP",3]), 3)
  dup.80 <- round(mean(samp.counts.80[samp.counts.80$GenGroup==i & samp.counts.80$Var2=="DUP",3]), 3)
  m <- cbind(i,del.dup, del.dup.80, del, del.80, dup, dup.80)
  means.gg <- rbind(means.gg, m)
}

means.gg <- as.data.frame(rbind(means.gg, cbind("All", round(mean(samp.counts$Freq)*2,3), 
            round(mean(samp.counts.80$Freq)*2,3),
            round(mean(samp.counts[samp.counts$Var2=="DUP", 3]),3),
            round(mean(samp.counts.80[samp.counts.80$Var2=="DUP", 3]),3),
            round(mean(samp.counts[samp.counts$Var2=="DEL", 3]),3),
            round(mean(samp.counts.80[samp.counts.80$Var2=="DEL", 3]),3))))

colnames(means.gg)[1] <- "GenGroup"

for (i in 2:ncol(means.gg)) {
  means.gg[,i] <- as.numeric(as.character(means.gg[,i]))
}

############ Perform a Mann-Whitney test ##########
kruskal.counts.del <- kruskal.test(x=samp.counts[samp.counts$Var2=="DEL",3], 
                                  g= as.factor(samp.counts[samp.counts$Var2=="DEL",11]))
capture.output(pairwise.wilcox.test(x=samp.counts[samp.counts$Var2=="DEL",3], 
                                  g=as.factor(samp.counts[samp.counts$Var2=="DEL",11]), 
                                  p.adjust.method="holm"), file = "../out/xhmm/counts-ind/counts_del_gg.posthoc")

kruskal.counts.dup <- kruskal.test(x=samp.counts[samp.counts$Var2=="DUP",3], 
                                   g= as.factor(samp.counts[samp.counts$Var2=="DUP",11]))
capture.output(pairwise.wilcox.test(x=samp.counts[samp.counts$Var2=="DUP",3], 
                                   g=as.factor(samp.counts[samp.counts$Var2=="DUP",11]), 
                                   p.adjust.method="holm"), file = "../out/xhmm/counts-ind/counts_dup_gg.posthoc")

########## Plots #############
ggplot(samp.counts, aes(Freq, fill = Var2)) + 
  geom_bar(alpha = 0.5, position = "identity") +
  xlab(label = "n CNVs per sample") +
  theme_bw(base_size = 29) 

ggplot(samp.counts, aes(Freq, fill = Var2)) + 
  geom_bar(alpha = 0.5, position = "identity") +
  geom_vline(data = means.gg[-9,], mapping=aes(xintercept = del), color = "#F8766D", size = 1.9) +
  geom_vline(data = means.gg[-9,], mapping=aes(xintercept = dup), color = "#00BFC4", size = 1.9) +
  xlab(label = "n CNVs per sample") + ylim(c(0,13.5)) +
  theme_bw(base_size = 29) +
  facet_grid(GenGroup~.)

ggplot(samp.counts.80, aes(Freq, fill = Var2)) + 
  geom_bar(alpha = 0.5, position = "identity") +
  geom_vline(data = means.gg[-9,], mapping=aes(xintercept = del.80), color = "#F8766D", size = 1.9) +
  geom_vline(data = means.gg[-9,], mapping=aes(xintercept = dup.80), color = "#00BFC4", size = 1.9) +
  xlab(label = "n CNVs per sample") +
  theme_bw(base_size = 29) + ylim(c(0,8)) +
  facet_grid(GenGroup~.)

##### Genetic groups origins
gg<- as.data.frame(table(samp.counts$Origin, samp.counts$GenGroup))
gg <- gg[gg$Freq!=0,]

##### Disecting 8 genetic group
## Canada, Mexico, Brazil and USA
eight <- droplevels(subset(samp.counts, GenGroup=="8"))
eight.countries <-droplevels(subset(samp.counts, Origin=="Canada" | Origin=="Mexico" | Origin=="Brazil" | Origin=="USA"))
gg8 <- as.data.frame(table(eight$Origin))
gg8$Freq2 <- gg8$Freq/2
gg8$Per <- (gg8$Freq2/sum(gg8$Freq2))*100
gg8.sh <- gg8[gg8$Freq2>2,c(1,3,4)]

ggplot(eight, aes(Freq, fill = Var2)) + 
  geom_bar(alpha = 0.5, position = "identity") +
  geom_vline(xintercept = mean(eight[eight$Var2=="DEL",3]), color = "#F8766D", size = 1.9) +
  geom_vline(xintercept = mean(eight[eight$Var2=="DUP",3]), color = "#00BFC4", size = 1.9) +
  xlab(label = "n CNVs per sample") + ylim(0,13) +
  theme_bw(base_size = 23)

ggplot(gg8.sh, aes(x="", y=Per, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  theme_bw(base_size = 27) + ylim(0,100) +
  coord_polar("y", start=0) +
  scale_fill_brewer(palette="Paired")
  #facet_wrap(~Var1)


