#library(xhmmScripts)
library(ggplot2)

#Load data
aux <- read.delim("../data/xhmm/LensXhmm.aux_xcnv")
cnv <- read.delim("../data/xhmm/LensXhmm.xcnv")
meta <- read.delim("../meta/LDP_meta.txt")

##Add metadata to cnv
cnv.meta <- numeric(0)
for (i in 1:nrow(meta)) {
  a <- cnv[cnv$SAMPLE == as.character(meta[i,1]),]
  b <- cbind(a, meta[i,-c(1,2)])
  cnv.meta <- rbind(cnv.meta, b)
}
####################### Depth ##############################
###### Mean depth according to sample
ind.depth <- cbind(aggregate(ORIG_RD ~ SAMPLE, aux, mean), aggregate(RD ~ SAMPLE,aux,  mean)[,2])

ggplot(ind.depth, aes(ORIG_RD)) + 
  geom_histogram(bins=50) 
ggplot(ind.depth, aes(ind.depth[,3])) + 
  geom_histogram(bins=50)

##### Depth site
#Original depth
aux$ori.rd.round <- round(aux$ORIG_RD)
aux$ori.rd.round[aux$ori.rd.round >= 200] <- 200

ggplot(aux, aes(ori.rd.round)) +
  geom_histogram(bins = 100) +
  theme_bw(base_size = 19)

#Normalized depth
ggplot(aux, aes(RD)) +
  geom_histogram(bins = 70) +
  theme_bw(base_size = 19)

#Amount of regions
variants <- droplevels(unique(cnv.ch[,c(2:6)]))
regions <- levels(variants$INTERVAL)

