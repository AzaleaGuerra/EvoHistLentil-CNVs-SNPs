#An evolutionary look into the history of lentil reveals unexpected diversity

This repository contains the script used in the article "An evolutionary look into the history of lentil reveals unexpected diversity"

##Prerequistes
XHMM
PLINK 1.07
VCFtools 0.1.15
Admixture 1.3
ADZE 1.0
Bayescan

###R packages
stringr
adegenet
SNPRelate
boa
coda
RIdeogram
PCAdapt
ggplot2

##Directories
**bin**: contains the scripts and is the working directory for them.
**meta**: information related to data collection, grouping, populations, etc.

##Data
Data is available at [https://knowpulse.usask.ca/AGILE/2](https://knowpulse.usask.ca/AGILE/2)  under Associated Datasets.

##Script organization
**CNV discovery**
CNV discovery was performed using XHMM (eXome-Hidden Markov Model). A detail explanation of this can be found in the file ``XHMM_lens.md``

**CNV**
Descriptors of the CNVs identified, including lenght, counts, private, and plotting: ``CNV-Basics.R``, ``CNV_Freq&priv.R``, ``CVN-sample-count.R``, ``CNV-idiogram_chr.R``, and ``CNV-length.R``.

Genes affected by CNVs: ``CNV-Genes.R`` and ``CNV-ReGs.R`` for resistance genes.


**Genetic diversity and differentiation**
Allele richness and private private alleles: ``adze.sh``
Vst (differentiation for CNV): ``Vst.R``
PCA of SNPs: ``PCA_LDP.R``
PCA using read depth of CNV: ``PCA-CNV-DP.R``


**Selection tests**
PCAdapt: ``PCAdapt.R``
Bayescan: ``Bayescan.sh`` and ``Bayescan_plot.sh``


##Contact
Azalea Guerrra Garcia
a.guerra@usask.ca
azalea.guerra@iecologia.unam.mx
azalea.guerra.g@gmail.com
