library(SNPRelate)
library(RColorBrewer)
library(ggplot2)
library(scatterplot3d)

###Tranform from vcf to gds
snpgdsVCF2GDS("../data/vcf/LDP_Lc2_MSP10_Q30_DP3_MAC5_CNV.vcf.gz", method = "copy.num.of.ref",
              out.fn="../data/gds/LDP_Lc2_MSP10_Q30_DP3_MAC5_CNV.gds")

#Load data
lentil <- snpgdsOpen("../data/gds/LDP_Lc2_MSP10_Q30_DP3_MAC10_NoCNV.gds")
sample.id <- read.gdsn(index.gdsn(lentil, "sample.id"))
meta <- read.delim("../meta/LDP_meta_SNPs.txt")

###PCA
pca <- snpgdsPCA(lentil, num.thread=3, autosome.only = FALSE, maf = 0.05)

# % of explined variation de variación contenido por los primeros componentes
pc.percent<- round(pca$varprop*100, 2)

## Order data according to genetic group in a data frame
pca.lentil <- data.frame(sample.id = pca$sample.id,
                     EV1 = pca$eigenvect[,1],  
                     EV2 = pca$eigenvect[,2],   
                     EV3 = pca$eigenvect[,3],
                     EV4 = pca$eigenvect[,4],
                     EV5 = pca$eigenvect[,5],
                     EV6 = pca$eigenvect[,6],
                     EV7 = pca$eigenvect[,7],
                     EV8 = pca$eigenvect[,8],
                     Pop = as.factor(meta$GenGroup),
                     Pop80 = as.factor(meta$GenGroup80),
                     stringsAsFactors = FALSE)

#Plot
my_palette <- brewer.pal(name="Set1", n=8)

ggplot(pca.lentil, aes(x=EV1, y=EV2)) +
  geom_point(aes(fill=Pop), size=10, shape=21, colour="black", alpha=0.8) + 
  xlab(paste0("Principal component 1 explains ", pc.percent[1], "%")) +
  ylab(paste0("Principal component 2 explains ", pc.percent[2], "%")) +
  scale_fill_manual(values=my_palette) + 
  theme_bw() + 
  theme(axis.text=element_text(size=30), axis.title=element_text(size=37), 
        legend.text=element_text(size=30), legend.title = element_blank(), 
        legend.key = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=7)))

palette(adjustcolor(brewer.pal(name="Set1", n=8), alpha.f = 0.8))
source('addgrids3d.r')

s3d <- scatterplot3d(x = pca.lentil$EV1,
              y = pca.lentil$EV3,
              z = pca.lentil$EV2,
              pch = "", cex.lab = 2.55, bg = pca.lentil$Pop, 
              cex.symbols = 1.5, cex.axis = 1.5,
              angle = 40, box=F, scale.y = 1,
              xlab = paste0("PC 1 explains ", pc.percent[1], "%"),
              ylab = paste0("PC 3 explains ", pc.percent[3], "%"),
              zlab = paste0("PC 2 explains ", pc.percent[2], "%"),
              xlim = c(-0.05,0.10), ylim = c(-0.05, 0.16), zlim = c(-0.10, 0.15))

addgrids3d(x = pca.lentil$EV1,
           y = pca.lentil$EV3,
           z = pca.lentil$EV2, grid = c("xy", "xz", "yz"),
           angle = 40, scale.y = 1,
           xlim = c(-0.1,0.10)) 

s3d$points3d(x = pca.lentil$EV1,
             y = pca.lentil$EV3,
             z = pca.lentil$EV2,
             pch = 21, bg = pca.lentil$Pop, cex = 3)

#######Violin plot
viol.data <- data.frame(rbind(cbind(pca.lentil$EV1, rep("PC1", nrow(pca.lentil)), as.character(pca.lentil$Pop), as.character(pca.lentil$Pop)), 
                              cbind(pca.lentil$EV2, rep("PC2", nrow(pca.lentil)), as.character(pca.lentil$Pop), as.character(pca.lentil$Pop)), 
                              cbind(pca.lentil$EV3, rep("PC3", nrow(pca.lentil)), as.character(pca.lentil$Pop), as.character(pca.lentil$Pop)),
                              cbind(pca.lentil$EV4, rep("PC4", nrow(pca.lentil)), as.character(pca.lentil$Pop), as.character(pca.lentil$Pop))))

viol.data$X1 <- as.numeric(as.character(viol.data$X1))

ggplot(viol.data, aes(x=X2, y=X1, fill = X4)) + 
  geom_violin(scale = "width", width = 0.7, trim = TRUE) +
  scale_fill_manual(values=my_palette) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), colour = "#88877f")+
  theme_bw() +
  theme(axis.text=element_text(size=18), axis.title= element_blank(), 
        legend.text=element_text(size=18), legend.title = element_blank(), 
        legend.key = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=7)))


viol.data80 <- data.frame(rbind(cbind(na.omit(pca.lentil)$EV1, rep("PC1", nrow(na.omit(pca.lentil))), as.character(na.omit(pca.lentil)$Pop), as.character(na.omit(pca.lentil)$Pop)), 
                              cbind(na.omit(pca.lentil)$EV2, rep("PC2", nrow(na.omit(pca.lentil))), as.character(na.omit(pca.lentil)$Pop), as.character(na.omit(pca.lentil)$Pop)), 
                              cbind(na.omit(pca.lentil)$EV3, rep("PC3", nrow(na.omit(pca.lentil))), as.character(na.omit(pca.lentil)$Pop), as.character(na.omit(pca.lentil)$Pop)),
                              cbind(na.omit(pca.lentil)$EV4, rep("PC4", nrow(na.omit(pca.lentil))), as.character(na.omit(pca.lentil)$Pop), as.character(na.omit(pca.lentil)$Pop))))

viol.data80$X1 <- as.numeric(as.character(viol.data80$X1))

ggplot(viol.data80, aes(x=X2, y=X1, fill = X4)) + 
  geom_violin(scale = "width", width = 0.7, trim = TRUE) +
  scale_fill_manual(values=my_palette) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), colour = "#88877f")
