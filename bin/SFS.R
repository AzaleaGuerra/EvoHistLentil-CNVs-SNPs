###Site Frequency Spectrum using allele count, considering all polymorphic SNPs per group
setwd("../data/freq-snps/")
library(ggplot2)

#Load list of file names
pop <- list.files(pattern = ".frq.count", full.names = FALSE)
#Load lifes
counts <- lapply(setNames(pop, make.names(gsub("*.frq.count$", "", pop))), read.delim)
setwd("../../bin/")
#Remove monomorphic sites within populations. 
counts <- lapply(counts,  function(x) { x <- x[x$C1!=0 & x$C2!=0,]; x})
#Add a column with pop info per data frame and a column with the least frequent allele count inside each pop
pop <- names(counts)

for (i in 1:length(counts)) {
  counts[[i]]$Pop <- rep(pop[i], nrow(counts[[i]]))
  counts[[i]]$C3 <- pmin(counts[[i]]$C1, counts[[i]]$C2)
}

#Dataframe containg the frequencies and counts
resum <- numeric(0)
for (i in 1:length(counts)) {
  for (j in 1:max(counts[[i]]$C2)) {
    aa <- cbind(pop[i], j, nrow(counts[[i]][counts[[i]]$C2 == j,]), nrow(counts[[i]][counts[[i]]$C2 == j,])/nrow(counts[[i]]))
    resum <- rbind(resum, aa)
  }
}

#Split list into dataframes
list2env(counts, envir = .GlobalEnv)

################# EXPECTED SFS ######################
#Get tetha values for each population
tetha <- numeric(0)

for (j in 1:length(counts)) {
  for (i in 1:(max(counts[[j]]$C2)-1)) {
    ind <- rbind(i, 1/i)
  }
  suma <- sum(ind)
  tetha <- rbind(tetha, cbind(pop[j], nrow(counts[[j]])/suma))
}
tetha <-as.data.frame(tetha)
tetha$V2 <- as.numeric(tetha$V2)

#Get SFS for each pop
ex.sfs <- list(0)
for (j in 1:length(counts)) {
  e.sfs <- numeric(0)
  for (i in 1:max(counts[[j]]$C2)) {
  b <- tetha[j,2]/i
  e.sfs <- rbind(e.sfs, cbind(i, b))
  }
  ex.sfs[[j]] <- as.data.frame(e.sfs)
  ex.sfs[[j]]$b <- as.numeric(as.character(ex.sfs[[j]]$b))
  ex.sfs[[j]]$i <- as.numeric(as.character(ex.sfs[[j]]$i))
  ex.sfs[[j]]$Pop <- pop[j]
  ex.sfs[[j]]$prop <- ex.sfs[[j]]$b/sum(ex.sfs[[j]]$b)
}
names(ex.sfs) <- paste(pop, "ex", sep = ".")
#lapply(ex.sfs, function(z) ggplot(z, aes(x=i, y = prop))+  geom_point() + geom_line()+theme_bw())
#Split list into dataframes
list2env(ex.sfs, envir = .GlobalEnv)

############### PLOT ##############
##Bind data according to populations status
lentil <- rbind(GG1, GG2, GG3, GG4, GG5, GG6, GG7, GG8)
lentil <- cbind(POS = lentil[,2], Type = rep("SNPs", nrow(lentil)), lentil[,c(1,7)], Count = lentil[,6])
lentil.ex <- rbind(GG1.ex, GG2.ex, GG3.ex, GG4.ex, GG5.ex, GG6.ex, GG7.ex, GG8.ex)
lentil.ex$Type <- rep("SNPs", nrow(lentil.ex)) 

ggplot(lentil, aes(x=Count, fill="", y = ..prop..)) + geom_bar(position = "identity") + 
  facet_grid(~Pop, scales = "free_x")+theme_bw(base_size = 27)+theme(legend.position="none") +
  scale_fill_manual(values = "#b8b7b7") +
  ylim(c(0,0.6)) + 
  geom_point(data = lentil.ex, aes(x=i, y= prop)) + geom_line(data = lentil.ex, aes(x=i, y= prop))
##Note about color of the bars: If a diffeent color per pop is needed, fill=Pop argument must be added in aes() argument

##Plot including the CNV
cnv <- read.delim("../data/CNVs/freq/cnv_freq2.txt", sep = "\t")
cnv$Type <- rep("CNV", nrow(cnv))

#Bind snps and cnvs
snps.cnv <- rbind(lentil, cnv[,c(1,2,3,6,7)])
snps.cnv$Type <- as.factor(snps.cnv$Type)

ggplot(snps.cnv[snps.cnv$Pop == "GG1" | snps.cnv$Pop == "GG2" | snps.cnv$Pop == "GG3" | snps.cnv$Pop == "GG4",] , 
       aes(x=Count, fill=Type, colour=Type, y = ..prop..)) + geom_bar(position = "identity", alpha=0.6) + 
  facet_grid(~Pop, scales = "free_x")+theme_bw(base_size = 27)+theme() +
  scale_fill_manual(values = c("#7dca75", "#8176cc", "#8176cc")) +
  scale_colour_manual(values = c("#7dca75", "#8176cc", "#8176cc")) +
  ylim(c(0,0.8)) +
  geom_point(data=lentil.ex[lentil.ex$Pop =="GG1" | lentil.ex$Pop =="GG2" | lentil.ex$Pop =="GG3" | lentil.ex$Pop =="GG4",],
             aes(x=i, y= prop), colour="black", size=1.1) + 
  geom_line(data=lentil.ex[lentil.ex$Pop=="GG1" | lentil.ex$Pop=="GG2" | lentil.ex$Pop=="GG3" | lentil.ex$Pop=="GG4",], 
            aes(x=i, y= prop), colour="black")


ggplot(snps.cnv[snps.cnv$Pop == "GG5" | snps.cnv$Pop == "GG6" | snps.cnv$Pop == "GG7" | snps.cnv$Pop == "GG8",] , 
       aes(x=Count, fill=Type, colour=Type, y = ..prop..)) + geom_bar(position = "identity", alpha=0.6) + 
  facet_grid(~Pop, scales = "free_x")+theme_bw(base_size = 27)+theme() +
  scale_fill_manual(values = c("#7dca75", "#8176cc", "#8176cc")) +
  scale_colour_manual(values = c("#7dca75", "#8176cc", "#8176cc")) +
  ylim(c(0,0.8)) 
  geom_point(data=lentil.ex[lentil.ex$Pop =="GG5" | lentil.ex$Pop =="GG6" | lentil.ex$Pop =="GG7" | lentil.ex$Pop =="GG8",],
             aes(x=i, y= prop), colour="black", size=1.1) + 
  geom_line(data=lentil.ex[lentil.ex$Pop=="GG5" | lentil.ex$Pop=="GG6" | lentil.ex$Pop=="GG7" | lentil.ex$Pop=="GG8",], 
            aes(x=i, y= prop), colour="black")

