library(ggplot2)
library(RColorBrewer)
library(plotly)
library(dplyr)
library(scatterplot3d)


dp <- read.delim("../data/xhmm/LensXhmm.filtered_centered.RD.txt", row.names = 1)
dp <- na.omit(dp)
dp <- dp[row.names(dp)!= "CN_105715-merged" & row.names(dp)!="QG-1-merged" & row.names(dp)!="SB-1-merged" &
           row.names(dp)!="964a-46" & row.names(dp)!="IG_72539" & row.names(dp)!="IG_72541" & row.names(dp)!="Lupa" & row.names(dp)!="PI_429838_LSP", ]

meta <- read.delim("../meta/LDP_meta.txt")

#keep only the depth of CNV regions present in more than 20% of the samples
dp.f <- numeric(0)
for (i in 1:ncol(dp)) {
  if (length(dp[dp[,i]>80 | dp[,i]< -80, i]) > 5) {
    dp.f <- cbind(dp.f, dp[,i])
  }
}
dp.f <- as.data.frame(dp.f)
row.names(dp.f) <- row.names(dp)
meanDP <- as.data.frame(rowMeans(dp.f))
row.names(meanDP) <- row.names(dp)

#PCA analysis
pca.dp <- prcomp(dp.f)
pcs <- as.data.frame(pca.dp$x)

#Bind meta data
pcs.meta <- numeric(0)
for (i in 1:nrow(meta)) {
  if(nrow(pcs[as.character(row.names(pcs)) == as.character(meta[i,1]),]) > 0) {
  a <- pcs[as.character(row.names(pcs)) == as.character(meta[i,1]),]
  b <- cbind(meta[i, c(1,8,9)], a)
  pcs.meta <- rbind(b, pcs.meta)
  }
}

#Estimate the explained variance
pca.var <- (pca.dp$sdev^2 / sum(pca.dp$sdev^2))*100
pca.var <- round(pca.var, digits = 2)

#Plot
my_palette <- brewer.pal(name="Set1", n=8)

ggplot(pcs.meta, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill=as.factor(GenGroup)), size=10, shape=21, colour="black", alpha=0.8) + 
  xlab(paste0("PC 1 explains ", pca.var[1], "%")) +
  ylab(paste0("PC 2 explains ", pca.var[2], "%")) +
  scale_fill_manual(values=my_palette) + 
  theme_bw() + 
  theme(axis.text=element_text(size=35), axis.title=element_text(size=37), 
        legend.text=element_text(size=35), legend.title = element_blank(), 
        legend.key = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=7)))

ggplot(pcs.meta, aes(x=PC2, y=PC3)) +
  geom_point(aes(fill=as.factor(GenGroup)), size=10, shape=21, colour="black", alpha=0.8) + 
  xlab(paste0("PC 2 explains ", pca.var[2], "%")) +
  ylab(paste0("PC 3 explains ", pca.var[3], "%")) +
  scale_fill_manual(values=my_palette) + 
  theme_bw() + 
  theme(axis.text=element_text(size=35), axis.title=element_text(size=37), 
        legend.text=element_text(size=35), legend.title = element_blank(), 
        legend.key = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=7)))

palette(adjustcolor(brewer.pal(name="Set1", n=8), alpha.f = 0.8))

source('addgrids3d.r')
s3d <- scatterplot3d(x = pcs.meta$PC2,
                     y = pcs.meta$PC3,
                     z = pcs.meta$PC4, pch = "", bg=pcs.meta$GenGroup,
                     grid=FALSE, box=FALSE,
                     cex.lab = 2, angle = 30, scale.y = 1.5, cex.symbols = 1.7, cex.axis = 1.2, 
                     xlab = "PC 1 explains 34%",
                     ylab = "PC 2 explains 10%",
                     zlab = "PC 3 explains 6%")
addgrids3d(x = pcs.meta$PC2,
           y = pcs.meta$PC3,
           z = pcs.meta$PC4, grid = c("xy", "xz", "yz"),
           angle = 30, scale.y = 1.5)
s3d$points3d(x = pcs.meta$PC2,
              y = pcs.meta$PC3,
              z = pcs.meta$PC4,
              pch = 21, bg = pcs.meta$GenGroup, cex = 1.7)
#legend("right", legend = levels(as.factor(pcs.meta$GenGroup)), pch = 21, bg = pcs.meta$GenGroup)

plot_ly(pcs.meta, x=~PC1, y=~PC2, z=~PC3, 
        color=~GenGroup, colors = my_palette) 

