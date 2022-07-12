#library(xhmmScripts)
library(ggplot2)
library(stringr)
library(RColorBrewer)

#Load data. Only CNVs within chromosomes
cnv<- read.delim("../data/xhmm/Lcul_CNV_chr.txt", sep = "\t")
meta <- read.delim("../meta/LDP_meta.txt")

#Keep sample with ancestryproportion > 0.80
cnv.80 <- droplevels(cnv[!is.na(cnv$GenGroup80),])

#######################################################################################
#### CNV Frequencies within genetic groups
cnv.freq <- as.data.frame(table(cnv$INTERVAL, cnv$CNV, cnv$GenGroup))
cnv.freq <- cnv.freq[cnv.freq$Freq != 0, ]
cnv.freq <- cbind(cnv.freq[c(1:2)], str_split_fixed(cnv.freq$Var1, pattern = ":", n=2), cnv.freq[c(3:4)])
cnv.freq <- cbind(cnv.freq[,c(1:3)], str_split_fixed(cnv.freq[,4], pattern = "-", n=2), cnv.freq[,c(5:6)])
colnames(cnv.freq) <- c("INTERVAL", "CNV", "CHR", "Start", "End",  "GenGroup", "Count")

#Frequencies 80% ancestry filter
cnv.freq80 <- as.data.frame(table(cnv$INTERVAL, cnv$CNV, cnv$GenGroup80))
cnv.freq80 <- cnv.freq80[cnv.freq80$Freq != 0, ]
cnv.freq80 <- cbind(cnv.freq80[c(1:2)], str_split_fixed(cnv.freq80$Var1, pattern = ":", n=2), cnv.freq80[c(3:4)])
cnv.freq80 <- cbind(cnv.freq80[,c(1:3)], str_split_fixed(cnv.freq80[,4], pattern = "-", n=2), cnv.freq80[,c(5:6)])
colnames(cnv.freq80) <- c("INTERVAL", "CNV", "CHR", "Start", "End",  "GenGroup80", "Count")

####################################################
################ PLOTS ######################

###### Plot DEL vs DUP
ggplot(cnv.freq) +
  geom_bar(aes(x=GenGroup, fill = CNV), position= "fill", size = 0.7, alpha=0.8) +
  theme_bw(base_size = 25)  

## CNV, >80%
ggplot(cnv.freq80) +
  geom_bar(aes(x=GenGroup80, fill = CNV), position= "fill", size = 0.7, alpha=0.8) +
  theme_bw(base_size = 25) 

##################################################
####### Plot frequencies
ggplot(cnv.freq, aes(Count, fill=CNV, y = ..prop.., colour=CNV)) +
  geom_bar(position = "identity", alpha=0.3) +
  facet_grid(~GenGroup, scales = "free_x") +
  theme_bw(base_size = 29) 

###############################################################
#### Identify private CNVs
priv <- row.names(as.matrix(which(table(cnv.freq$INTERVAL)==1)))

#Private CNVs, 80% ancestry filter
priv80 <- row.names(as.matrix(which(table(cnv.freq80$INTERVAL)==1)))

##################################################
####### Plot frequencies
ggplot(cnv.freq, aes(Count, fill=CNV, y = ..prop.., colour=CNV)) +
  geom_bar(position = "identity", alpha=0.3) +
  facet_grid(~GenGroup, scales = "free_x") +
  theme_bw(base_size = 29) 

cnv.f <- as.data.frame(table(cnv$INTERVAL, cnv$GenGroup))
cnv.f <- cnv.f[cnv.f$Freq != 0, ]
cnv.f <- cbind(cnv.f, str_split_fixed(cnv.f$Var1, pattern = ":", n=2), paste(cnv.f$Var1, cnv.f$Var2, sep = "-"))
colnames(cnv.f) <- c("INTERVAL", "GenGroup", "Freq", "CHR", "Pos", "Inter.CNV")

ggplot(cnv.f, aes(Freq, y = ..prop..)) +
  geom_bar(position = "identity", alpha=0.7, fill="#777acd", color="#777acd") +
  facet_grid(~GenGroup, scales = "free_x") +
  theme_bw(base_size = 29) 
