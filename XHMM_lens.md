#Structral Variant Discovery

Scripts and files related to this process are in ``CNVdiscovery`` directory:
``CNVdiscovery``:
--``bin``: scripts and executable files are located here. This is the **working directory**.
--``bin/nohup_out``: logs and error messages files
--``data``
--``meta``: samples names
--``out``: contains xhmm outputs

##XHMM
 [XHMM](https://www.sciencedirect.com/science/article/pii/S000292971200417X?via%3Dihub): Caller focused on exome data. It uses PCA to normalize exome read depth and a hidden Markov model (HMM) to discover exon-resolution CNV and genotype variation across samples.
 [XHMM tutorial](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4065038/)

**Steps:**
**1** Prepare data
**1.1** Add Group parameter to bam files (GATK requirement). This was performed with Picard 2.22.1. **Done**.
Script `AddGroupRead.sh`:
```
#!/bin/bash
 
bin="./"
bam="../data/LDP-bam"
data="../data/DepthCov"
meta="../meta"
ref="../data/Lculinaris/Lens_culinaris_2.0.fasta"

while read prefix
do
 if [ ! -f $bam/$prefix ]
then

java -Xmx100g -jar picard.jar AddOrReplaceReadGroups \
      I=$bam/$prefix.best-sorted-rmdup.bam \
      O=$bam/bam-sort/$prefix.sort.bam \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=lane \
      RGSM=$prefix

fi
done < ../meta/IDs.txt
```

**1.2**  Merge bam files with sambamba 0.7.1. **Done**. 
Fragment of the script `merge-bams.sh`:
```
#!/bin/bash

ncores="20"
bin="./"
bam="../data/LDP-bam/bam-sort"
meta="../meta"
files=echo

./sambamba/bin/sambamba-0.7.1 merge	$bam/bams_merged/1-30.bam $bam/CN_105715-merged.sort.bam \
$bam/ILL_4609.sort.bam \
$bam/IG_1959.sort.bam \
$bam/CDC_GREENSTAR.sort.bam \
$bam/SHASTA.sort.bam \
$bam/CDC_KR-1.sort.bam \
$bam/PI_426797_LSP.sort.bam \
$bam/ILL_8007.sort.bam \
$bam/ILL_9.sort.bam \
$bam/IG_1706.sort.bam \
$bam/CDC_MAXIM.sort.bam \
$bam/ILL_213.sort.bam \
$bam/ILL_5722.sort.bam \
$bam/PI_431884_LSP.sort.bam \
$bam/CN_105789.sort.bam \
$bam/PI_432245_LSP.sort.bam \
$bam/PI_182217_LSP.sort.bam \
$bam/PI_432005_LSP.sort.bam \
$bam/ILL_618.sort.bam \
$bam/ILL_1983.sort.bam \
$bam/CN_105895.sort.bam \
$bam/ILL_2507.sort.bam \
$bam/ILL_3347.sort.bam \
$bam/ILL_3597.sort.bam \
$bam/ILL_4768.sort.bam \
$bam/ILL_7946.sort.bam \
$bam/ILL_975.sort.bam \
$bam/ILL_5151.sort.bam \
$bam/ILL_5480.sort.bam \
$bam/ILL_5883.sort.bam \
	-t $ncores
```

**1.3** Create a CSI index. FAI index do not support chromosomes longer than 530 Mbp and lentil Chr2 is ~610 Mbp. Script ``index.sh``.
**Done**.

**1.4** Create variant dictionary (Picard 2.22.1) and index (Samtools 1.9),  script ``VariantDictionary.sh`` **Done**.
```
#!/bin/bash
 
bin="./"
bam="../data/LDP-bam/bam-sort"
data="../data"
meta="../meta"
ref="../data/Lculinaris/Lens_culinaris_2.0.fasta"

#Create a variant dictonary using Picard
java -jar picard.jar CreateSequenceDictionary R=$ref O=../data/Lculinaris/Lens_culinaris_2.0.dict

#Create index with samtools
samtools faidx $ref
```

**1.5** Estimate depth of coverag with GATK 4.1.8.1; scripts ``DepthOfCoverage4_1-90``, ``DepthOfCoverage4_91-210`` and ``DepthOfCoverage4_211-329``. 
In order to run in parallel this step, which is the most time consuming,   samples (bam files) were splitted. 
  
```
bin="./"
bam="../data/LDP-bam/bam-sort/bams_merged"
data="../data"
meta="../meta"
gatk="$bin/gatk-4.1.8.1"
ref="../data/Lculinaris/Lens_culinaris_2.0.fasta"
genes="../data/Lculinaris/Lens_culinaris_2.0.genes.bed"

for prefix in 1-30 31-60 61-90; do

./$gatk/gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx55g" \
	DepthOfCoverage \
	-R $ref \
    -L $genes\
    -I $bam/$prefix.bam \
    -O $data/DepthCov/$prefix.4data \
    --omit-depth-output-at-each-base \
    --omit-locus-table \
    --include-ref-n-sites \
    --include-deletions ;

done 
```

GATK 4 uses "," to separte columns. Replace "," for tab. 
```
for prefix in 1-30 31-60 61-90 91-120 121-150 151-180 181-210 211-240 241-270 271-300 301-329; do
	sed 's/,/\t/g' ../data/DepthCov/${prefix}.4data.sample_interval_summary > ../data/DepthCov/${prefix}.data.sample_interval_summary
done
```


**2** Identification and filtering of CNV with XHMM 1.0.
**2.1** `xhmm-MergeDepth.sh` **Done**
```
#!/bin/bash
 
ncores="30"
bin="./"
data="../data"
meta="../meta"
out="../out"

#Combines GATK Depth-of-Coverage outputs for multiple samples (at same loci):
./xhmm/xhmm --mergeGATKdepths -o $out/xhmm/LensXhmm.RD.txt \
--GATKdepths $data/DepthCov/1-30.data.sample_interval_summary \
--GATKdepths $data/DepthCov/31-60.data.sample_interval_summary \
--GATKdepths $data/DepthCov/61-90.data.sample_interval_summary \
--GATKdepths $data/DepthCov/91-120.data.sample_interval_summary \
--GATKdepths $data/DepthCov/121-150.data.sample_interval_summary \
--GATKdepths $data/DepthCov/151-180.data.sample_interval_summary \
--GATKdepths $data/DepthCov/181-210.data.sample_interval_summary \
--GATKdepths $data/DepthCov/211-240.data.sample_interval_summary \
--GATKdepths $data/DepthCov/241-270.data.sample_interval_summary \
--GATKdepths $data/DepthCov/271-300.data.sample_interval_summary \
--GATKdepths $data/DepthCov/301-329.data.sample_interval_summary

```
**2.2** `xhmm-filtersPCA.sh` **Done**
```
#!/bin/bash
 
ncores="30"
bin="./"
data="../data"
meta="../meta"
out="../out"

#Filters samples and targets and then mean-centers the targets:
./xhmm/xhmm --matrix -r $out/xhmm/LensXhmm.RD.txt \
--centerData --centerType target \
-o $out/xhmm/LensXhmm.filtered_centered.RD.txt \
--outputExcludedTargets $out/xhmm/LensXhmm.filtered_centered.RD.txt.filtered_targets.txt \
--outputExcludedSamples $out/xhmm/LensXhmm.filtered_centered.RD.txt.filtered_samples.txt \
--minTargetSize 10 --maxMeanTargetRD 3000 \
--minMeanTargetRD 3

#Runs PCA on mean-centered data:
./xhmm/xhmm --PCA -r $out/xhmm/LensXhmm.filtered_centered.RD.txt --PCAfiles $out/xhmm/LensXhmm.RD_PCA

#Normalizes mean-centered data using PCA information:
./xhmm/xhmm --normalize -r $out/xhmm/LensXhmm.filtered_centered.RD.txt --PCAfiles $out/xhmm/LensXhmm.RD_PCA \
--normalizeOutput $out/xhmm/LensXhmm.PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

#Filters and z-score centers (by sample) the PCA-normalized data:
./xhmm/xhmm --matrix -r $out/xhmm/LensXhmm.PCA_normalized.txt \
--centerData --centerType sample --zScoreData \
-o $out/xhmm/LensXhmm.PCA_normalized.filtered.sample_zscores.RD.txt \
--outputExcludedTargets $out/xhmm/LensXhmm.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--outputExcludedSamples $out/xhmm/LensXhmm.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt 
#--maxSdTargetRD 30

#Filters original read-depth data to be the same as filtered, normalized data:
./xhmm/xhmm --matrix -r $out/xhmm/LensXhmm.RD.txt \
--excludeTargets $out/xhmm/LensXhmm.filtered_centered.RD.txt.filtered_targets.txt \
--excludeTargets $out/xhmm/LensXhmm.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--excludeSamples $out/xhmm/LensXhmm.filtered_centered.RD.txt.filtered_samples.txt \
--excludeSamples $out/xhmm/LensXhmm.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
-o $out/xhmm/LensXhmm.same_filtered.RD.txt

```

**2.3** `xhmm-discoverCNV.sh` **Done**
```
#!/bin/bash
 
ncores="30"
bin="./"
data="../data"
meta="../meta"
out="../out"
ref="../data/Lculinaris/Lens_culinaris_2.0.fasta"

#Discovers CNVs in normalized data:
./xhmm/xhmm --discover -p ./xhmm/params.txt \
-r $out/xhmm/LensXhmm.PCA_normalized.filtered.sample_zscores.RD.txt \
-R $out/xhmm/LensXhmm.same_filtered.RD.txt \
-c $out/xhmm/LensXhmm.xcnv -a $out/xhmm/LensXhmm.aux_xcnv

#Genotypes discovered CNVs in all samples:
./xhmm/xhmm --genotype -p ./xhmm/params.txt \
-r $out/xhmm/LensXhmm.PCA_normalized.filtered.sample_zscores.RD.txt \
-R $out/xhmm/LensXhmm.same_filtered.RD.txt \
-g $out/xhmm/LensXhmm.xcnv \
-F $ref \
-v $out/xhmm/LensXhmm.vcf
```

