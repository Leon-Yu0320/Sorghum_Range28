#  **Range 28 Sorghum epitranscriptome project**\

**Author: Li'ang Yu, Andrew Nelson**\
**Date: June 14th, 2022**

<img width="936" alt="Screenshot 2025-02-01 at 8 22 46 PM" src="https://github.com/user-attachments/assets/31bdbe6d-d5b6-4cca-af4e-9b75ece7bddd" />

<img width="1031" alt="Screenshot 2025-02-01 at 8 23 21 PM" src="https://github.com/user-attachments/assets/b03b3cc8-55a3-44f1-8d09-13c92818b8c9" />

- [**Range 28 Sorghum epitranscriptome project**](#range-28-sorghum-epitranscriptome-project)
- [**Section I Genomic variants and pan-genome analysis of sorghum accessions**](#section-i-genomic-variants-and-pan-genome-analysis-of-sorghum-accessions)
  - [**Download sorghum association pannel variants data**](#download-sorghum-association-pannel-variants-data)
    - [install ascp packages](#install-ascp-packages)
    - [Download the re-sequencing WGS reads](#download-the-re-sequencing-wgs-reads)
    - [Download sorghum variants calling file (VCF) based on selected panels](#download-sorghum-variants-calling-file-vcf-based-on-selected-panels)
  - [**Perform the snpEffect annotation**](#perform-the-snpeffect-annotation)
    - [Prepare the database](#prepare-the-database)
    - [Start with the annotation using database constructed above](#start-with-the-annotation-using-database-constructed-above)
  - [**Pan-genome consturction of the six accessions**](#pan-genome-consturction-of-the-six-accessions)
    - [Merge sub-reads for six accessions](#merge-sub-reads-for-six-accessions)
    - [Trim adaptor sequences (Treat as paired end)](#trim-adaptor-sequences-treat-as-paired-end)
    - [Genome asembly using Megahit pipeline](#genome-asembly-using-megahit-pipeline)
    - [Evaluation of genome assembly by BUSCO](#evaluation-of-genome-assembly-by-busco)
    - [Evaluation of genome assembly by QUAST](#evaluation-of-genome-assembly-by-quast)
    - [Extract unaligned sequences](#extract-unaligned-sequences)
    - [Remove duplicates](#remove-duplicates)
    - [Remove contamination](#remove-contamination)
      - [Download the NCBI nr/nt database](#download-the-ncbi-nrnt-database)
      - [Perform the cleaning steps](#perform-the-cleaning-steps)
    - [Novel gene annotation and prediciton](#novel-gene-annotation-and-prediciton)
      - [Using lift-off mapping strategy](#using-lift-off-mapping-strategy)
      - [Using RNA-seq reads mapping strategy](#using-rna-seq-reads-mapping-strategy)
      - [Mapping reads to Clean Pan sequences](#mapping-reads-to-clean-pan-sequences)
- [**Section II Gene expression and regulatory network analysis of experiments**](#section-ii-gene-expression-and-regulatory-network-analysis-of-experiments)
  - [STEP 1 Quantification of reads counts using featureCounts](#step-1-quantification-of-reads-counts-using-featurecounts)
    - [**Extraction of the primary transcripts and longest transcripts for lncRNAs based on sorghum gene annotation**](#extraction-of-the-primary-transcripts-and-longest-transcripts-for-lncrnas-based-on-sorghum-gene-annotation)
    - [**Perform the featureCounts for coding genes at exon level**](#perform-the-featurecounts-for-coding-genes-at-exon-level)
    - [check if meta-feature and features counted correctly for transcripts](#check-if-meta-feature-and-features-counted-correctly-for-transcripts)
    - [**Perform the featureCounts for lncRNAs (exon level)**](#perform-the-featurecounts-for-lncrnas-exon-level)
    - [check if meta-feature and features counted correctly for lncRNAs](#check-if-meta-feature-and-features-counted-correctly-for-lncrnas)
    - [**Perform the featureCounts for modififed sites (per PTM sites level)**](#perform-the-featurecounts-for-modififed-sites-per-ptm-sites-level)
    - [check if meta-feature and features counted correctly for PTMs](#check-if-meta-feature-and-features-counted-correctly-for-ptms)
  - [STEP 2 Normalization of read counts and differentially expressed genes (DEGs) analysis](#step-2-normalization-of-read-counts-and-differentially-expressed-genes-degs-analysis)
    - [**Install or load packages**](#install-or-load-packages)
    - [**Import FeatureCounts data**](#import-featurecounts-data)
    - [**Experimental design and DEGs analysis**](#experimental-design-and-degs-analysis)
    - [**Export data for inter-accession comparison**](#export-data-for-inter-accession-comparison)
    - [**Extraction of normalized counts**](#extraction-of-normalized-counts)
    - [**PCA based on vst transfomration**](#pca-based-on-vst-transfomration)
    - [**Distance Matrix clustering under WL and WW condition**](#distance-matrix-clustering-under-wl-and-ww-condition)
    - [**Transform into TPM**](#transform-into-tpm)
    - [**Pearson correaltion for each two samples based on normalized featureCounts/TPM**](#pearson-correaltion-for-each-two-samples-based-on-normalized-featurecountstpm)
  - [STEP 4 WGCNA and regulatory network analysis](#step-4-wgcna-and-regulatory-network-analysis)
    - [**Optimization of the nework construction**](#optimization-of-the-nework-construction)
      - [Prepare the PARAMETER_TABLE](#prepare-the-parameter_table)
      - [Prepare the R-code](#prepare-the-r-code)
      - [Prepare the bash file for loop Rcode](#prepare-the-bash-file-for-loop-rcode)
    - [**Prepare the Formal network construction**](#prepare-the-formal-network-construction)
      - [Load requried libraries and data input](#load-requried-libraries-and-data-input)
      - [Load expression data and meta-table](#load-expression-data-and-meta-table)
      - [Load phenotypic and metablomic data](#load-phenotypic-and-metablomic-data)
      - [Load genes with high MAD score](#load-genes-with-high-mad-score)
      - [Calculate the soft thresholding power (beta)](#calculate-the-soft-thresholding-power-beta)
      - [Step-by-step network construction and module detection](#step-by-step-network-construction-and-module-detection)
      - [Merging of modules whose expression profiles are very similar](#merging-of-modules-whose-expression-profiles-are-very-similar)
      - [Compare the epigengene values across samples from each module](#compare-the-epigengene-values-across-samples-from-each-module)
      - [Module-trait membership analysis](#module-trait-membership-analysis)
      - [Oxdative stress data correlation with expression](#oxdative-stress-data-correlation-with-expression)
      - [Build correlation between LI6800 data and expression module (morning data/afternoon data)](#build-correlation-between-li6800-data-and-expression-module-morning-dataafternoon-data)
      - [Calculate the gene-module membership (MM)](#calculate-the-gene-module-membership-mm)
      - [GS_MM based on particular trait](#gs_mm-based-on-particular-trait)
      - [Export network for cytoscape](#export-network-for-cytoscape)
    - [**Slection of trait-associated genes**](#slection-of-trait-associated-genes)
      - [Load MM and GS files](#load-mm-and-gs-files)
      - [Load gene annotation information](#load-gene-annotation-information)
      - [Load metatable containing the interested module-trait combinations](#load-metatable-containing-the-interested-module-trait-combinations)
      - [Plot MM-GS scater plot with cutoff marked](#plot-mm-gs-scater-plot-with-cutoff-marked)
    - [**Filter the weight-score of edge data and extract nodes**](#filter-the-weight-score-of-edge-data-and-extract-nodes)
      - [Match module membership (MM) to filtered nodes](#match-module-membership-mm-to-filtered-nodes)
  - [STEP 5 GO and KEGG enrichment analysis for different sets of genes](#step-5-go-and-kegg-enrichment-analysis-for-different-sets-of-genes)
    - [**Install and load packages**](#install-and-load-packages)
    - [**Load database**](#load-database)
    - [**Load different set of genes**](#load-different-set-of-genes)
    - [**Wrap gene function information to each module**](#wrap-gene-function-information-to-each-module)
    - [Perform functional enrichment](#perform-functional-enrichment)
    - [**Load KEGG pathway data**](#load-kegg-pathway-data)
    - [**Perform the KEGG enrichement analysis**](#perform-the-kegg-enrichement-analysis)
  - [STEP 6 Enrichment test of binding motifs from promoter regions](#step-6-enrichment-test-of-binding-motifs-from-promoter-regions)
    - [**Directories information for STEP 6**](#directories-information-for-step-6)
    - [**Extract promoter sequences baesd on sorghum annotation**](#extract-promoter-sequences-baesd-on-sorghum-annotation)
    - [**Replace the transcript ID to coordinate list**](#replace-the-transcript-id-to-coordinate-list)
    - [**Perform the enrichment analyiss based on AME pipeline**](#perform-the-enrichment-analyiss-based-on-ame-pipeline)
- [**Section III High-throughput data processing for post-transcriptional modifications (PTMs) using RNA-seq data**](#section-iii-high-throughput-data-processing-for-post-transcriptional-modifications-ptms-using-rna-seq-data)
  - [STEP 1 Identificaiton of PTMs using Tophat2 mapping and HAMR](#step-1-identificaiton-of-ptms-using-tophat2-mapping-and-hamr)
    - [**Container information**](#container-information)
    - [**Directories information for STEP 1**](#directories-information-for-step-1)
    - [**Index the reference genome**](#index-the-reference-genome)
    - [**Use trimmonmatic to remove adapters and low-quality reads**](#use-trimmonmatic-to-remove-adapters-and-low-quality-reads)
    - [**Check the LIB types before alignment**](#check-the-lib-types-before-alignment)
    - [**Perform reads alignment by Tophat2**](#perform-reads-alignment-by-tophat2)
    - [**Modify the format of bam files**](#modify-the-format-of-bam-files)
    - [**Perform PTMs detection by HAMR**](#perform-ptms-detection-by-hamr)
    - [**Perform PTMs detection by MODtect**](#perform-ptms-detection-by-modtect)
    - [**Perform the intersection among two replicates of two methods**](#perform-the-intersection-among-two-replicates-of-two-methods)
    - [**Annotate genes associated with HAMR modified nucleotides**](#annotate-genes-associated-with-hamr-modified-nucleotides)
    - [**Annotate genes associated with Modtect modified nucleotides**](#annotate-genes-associated-with-modtect-modified-nucleotides)
    - [**Intersetion of resutls of modified sites from HAMR and Modtect**](#intersetion-of-resutls-of-modified-sites-from-hamr-and-modtect)
    - [**Annotate genes associated with intersected modified nucleotides**](#annotate-genes-associated-with-intersected-modified-nucleotides)
    - [Subtract the information of mapping depth from HAMR](#subtract-the-information-of-mapping-depth-from-hamr)
    - [Subtract the information of mapping depth](#subtract-the-information-of-mapping-depth)
    - [**Directories information for STEP 2**](#directories-information-for-step-2)
    - [**Perform assignement of modified sites (PTMs) to respective adjacent primary transcripts**](#perform-assignement-of-modified-sites-ptms-to-respective-adjacent-primary-transcripts)
  - [STEP 3 Characterization of distribution of PTMs across diffrent genomic context and plot the meta-plot](#step-3-characterization-of-distribution-of-ptms-across-diffrent-genomic-context-and-plot-the-meta-plot)
    - [**Directories information for STEP 3**](#directories-information-for-step-3)
    - [**Generate a bed file used for parsing genomic-element feature (gff2bed.R)**](#generate-a-bed-file-used-for-parsing-genomic-element-feature-gff2bedr)
    - [**Genrate metafeature count for each sample**](#genrate-metafeature-count-for-each-sample)
    - [**Plot metaPlot of each sample or merged samples**](#plot-metaplot-of-each-sample-or-merged-samples)
  - [STEP 4 Qunatificaiton of modified genes and modified sites across transcriptome](#step-4-qunatificaiton-of-modified-genes-and-modified-sites-across-transcriptome)
    - [**Reformat the structure of PTMs count**](#reformat-the-structure-of-ptms-count)
    - [**Classify total numbers of mods and total genes been modifed**](#classify-total-numbers-of-mods-and-total-genes-been-modifed)
    - [**Plot overall patterns of modifed sites across samples**](#plot-overall-patterns-of-modifed-sites-across-samples)
    - [**Characterize the genes been modified genes from two relicas**](#characterize-the-genes-been-modified-genes-from-two-relicas)
    - [**Export numbers of PTMs  per gene based on averaging two reps**](#export-numbers-of-ptms--per-gene-based-on-averaging-two-reps)

# **Section I Genomic variants and pan-genome analysis of sorghum accessions**
## **Download sorghum association pannel variants data**
### install ascp packages
The installation and [**configuration**](https://www.jianshu.com/p/df34b99cbd43) can be accessed
```bash
cd /home/liangyu/3_software

###download the packages
wget http://download.asperasoft.com/download/sw/connect/3.7.4/aspera-connect-3.7.4.147727-linux-64.tar.gz

###Install the package
bash aspera-connect-3.7.4.147727-linux-64.sh

###add package into the PTH
echo 'export PATH=/home/liangyu/.aspera/connect/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

### Download the re-sequencing WGS reads
sample list from ENA database (examples)
```
/vol1/fastq/ERR809/002/ERR8097262/ERR8097262_1.fastq.gz /vol1/fastq/ERR809/002/ERR8097262/ERR8097262_2.fastq.gz
/vol1/fastq/ERR809/002/ERR8097502/ERR8097502_1.fastq.gz	/vol1/fastq/ERR809/002/ERR8097502/ERR8097502_2.fastq.gz
/vol1/fastq/ERR809/008/ERR8097508/ERR8097508_1.fastq.gz	/vol1/fastq/ERR809/008/ERR8097508/ERR8097508_2.fastq.gz
/vol1/fastq/ERR809/003/ERR8097613/ERR8097613_1.fastq.gz	/vol1/fastq/ERR809/003/ERR8097613/ERR8097613_2.fastq.gz
/vol1/fastq/ERR809/001/ERR8097641/ERR8097641_1.fastq.gz	/vol1/fastq/ERR809/001/ERR8097641/ERR8097641_2.fastq.gz
/vol1/fastq/ERR810/009/ERR8106149/ERR8106149_1.fastq.gz	/vol1/fastq/ERR810/009/ERR8106149/ERR8106149_2.fastq.gz
/vol1/fastq/ERR810/004/ERR8106214/ERR8106214_1.fastq.gz	/vol1/fastq/ERR810/004/ERR8106214/ERR8106214_2.fastq.gz
/vol1/fastq/ERR801/008/ERR8017558/ERR8017558_1.fastq.gz	/vol1/fastq/ERR801/008/ERR8017558/ERR8017558_2.fastq.gz
/vol1/fastq/ERR805/006/ERR8056536/ERR8056536_1.fastq.gz	/vol1/fastq/ERR805/006/ERR8056536/ERR8056536_2.fastq.gz
/vol1/fastq/ERR805/003/ERR8058553/ERR8058553_1.fastq.gz	/vol1/fastq/ERR805/003/ERR8058553/ERR8058553_2.fastq.gz
/vol1/fastq/ERR806/009/ERR8060479/ERR8060479_1.fastq.gz	/vol1/fastq/ERR806/009/ERR8060479/ERR8060479_2.fastq.gz
/vol1/fastq/ERR806/000/ERR8060670/ERR8060670_1.fastq.gz	/vol1/fastq/ERR806/000/ERR8060670/ERR8060670_2.fastq.gz
/vol1/fastq/ERR809/009/ERR8093059/ERR8093059_1.fastq.gz	/vol1/fastq/ERR809/009/ERR8093059/ERR8093059_2.fastq.gz
```
```bash
work_dir="/mnt/Knives/1_data/Sorghum_Pangenome"
IFS=$'\n';

for LINE in $(cat $work_dir/sample.list);
do
	pair1=$(echo ${LINE} | awk '{ print $1}')
	pair2=$(echo ${LINE} | awk '{ print $2}')

    ### forward reads download
	ascp -k 1 -QT -l 500m -P33001 \
    	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    	era-fasp@fasp.sra.ebi.ac.uk:$pair1\
    	/mnt/Knives/1_data/Sorghum_Pangenome/

    ### reverse reads download
	ascp -k 1 -QT -l 500m -P33001 \
    	-i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    	era-fasp@fasp.sra.ebi.ac.uk:$pair2\
    	/mnt/Knives/1_data/Sorghum_Pangenome/
done

### Count reads numbers for summary
for i in $(ls $work_dir/*.fastq.gz);
do
	echo $(zcat $i|wc -l)/4|bc >> read.counts
done
```

### Download sorghum variants calling file (VCF) based on selected panels 
```bash
genome_dir="/mnt/Milly/Work/Sorghum_epitranscriptomics_project/phytozome/Sbicolor/v3.1.1/assembly"
snpeff_dir="/home/liangyu/3_software/snpEff"
work_dir="/mnt/Knives/1_data/Sorghum_SNP/"

ascp -k 1 -QT -l 500m -P33001 \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    era-fasp@fasp.sra.ebi.ac.uk:/vol1/ERZ737/ERZ7374343/SAP.vcf.gz \
    /mnt/Knives/1_data/Sorghum_SNP/

###unzip samples
gunzip /mnt/Knives/1_data/Sorghum_SNP/SAP.vcf.gz

###Count numbers of variants across 400 samples (43,811,788 variants)
grep -v "##" /mnt/Knives/1_data/Sorghum_SNP/SAP.vcf | wc -l
```
**The accession panel list are shown as below:**
```
PI533871
PI533961
PI656041
PI656053
PI656076
PI656096
```
```bash
#extract depth (DP) for six samples from panel
vcftools --vcf /mnt/Knives/1_data/Sorghum_SNP/SAP.vcf \
    --keep /mnt/Knives/1_data/Sorghum_SNP/panel.list \
    --site-mean-depth && \
    mv /mnt/Knives/1_data/Sorghum_SNP/out.ldepth.mean \
    /mnt/Knives/1_data/Sorghum_SNP/panel_ldegpth.mean 

#suset the entire vcf information for six samples from panel (support multiple threads)
bcftools view -S /mnt/Knives/1_data/Sorghum_SNP/panel.list \
     /mnt/Knives/1_data/Sorghum_SNP/SAP.vcf \
    -o /mnt/Knives/1_data/Sorghum_SNP/select.vcf

#remove sites without polymorphysims (kept 8725449 out of a possible 43699819 Sites)
vcftools --vcf /mnt/Knives/1_data/Sorghum_SNP/select.vcf \
    --mac 1 \
    --recode --recode-INFO-all \
    --out filter_mac1_select
```

## **Perform the snpEffect annotation**
### Prepare the database 
Evaluate the SNP effects using [**SnpEff**](http://pcingola.github.io/SnpEff/se_running/) \
1: Build the [databases](https://www.cnblogs.com/zhanmaomao/p/10964636.html)

```bash
mkdir $snpeff_dir/data
mkdir $snpeff_dir/data/Sbicolor  
mkdir $snpeff_dir/data/genomes

###prepare the reference genome
ln -s $genome_dir/Sbicolor.fa  $snpeff_dir/data/genomes

###Prepare the gff files 
ln -s /home/liangyu/1_Project/3_RNA-mods/8_gffcompare/Sorghum.gff_primary_transcripts.gff3 $snpeff_dir/data/Sbicolor

###Build library 
echo "Sbicolor.genome : Sbicolor" >> $snpeff_dir/snpEff.config 

###Perform database constructions
java -Xmx16G -jar /home/liangyu/3_software/snpEff/snpEff.jar \
    build -gff3 Sbicolor

cd $snpeff_dir
### use cds and gff to create database
java -Xmx64g -jar snpEff.jar build -gff3 -v Sbicolor -noCheckProtein

###Check the database build from Sbicolor file
ls $snpeff_dir/data/Sbicolor

###Reformat sequences from small contigs to sequences for further looping steps
mv $snpeff_dir/sequence.bin  $snpeff_dir/sequence0.bin
```
Results should be like this:
```
drwxrwxrwx 4 liangyu liangyu      4096 Jul  8 16:41 ../
-rw-r--r-- 1 liangyu liangyu  40855612 Jul 11 11:49 cds.fa
-rw-r--r-- 1 liangyu liangyu  33330243 Jul  8 17:18 genes.gff
-rw-r--r-- 1 liangyu liangyu       643 Jul 11 11:46 protein.fa
-rw-r--r-- 1 liangyu liangyu   6170805 Jul 11 12:02 sequence.1.bin
-rw-r--r-- 1 liangyu liangyu   3107552 Jul 11 12:02 sequence.10.bin
-rw-r--r-- 1 liangyu liangyu   4677923 Jul 11 12:02 sequence.2.bin
-rw-r--r-- 1 liangyu liangyu   4947833 Jul 11 12:02 sequence.3.bin
-rw-r--r-- 1 liangyu liangyu   4109711 Jul 11 12:02 sequence.4.bin
-rw-r--r-- 1 liangyu liangyu   2435372 Jul 11 12:02 sequence.5.bin
-rw-r--r-- 1 liangyu liangyu   3184420 Jul 11 12:02 sequence.6.bin
-rw-r--r-- 1 liangyu liangyu   2490379 Jul 11 12:02 sequence.7.bin
-rw-r--r-- 1 liangyu liangyu   2217577 Jul 11 12:03 sequence.8.bin
-rw-r--r-- 1 liangyu liangyu   2866726 Jul 11 12:03 sequence.9.bin
-rw-r--r-- 1 liangyu liangyu     73556 Jul 11 12:03 sequence0.bin
-rwxr-xr-x 1 liangyu liangyu 720688555 Jul  8 17:18 sequences.fa*
-rw-r--r-- 1 liangyu liangyu  18313758 Jul 11 12:02 snpEffectPredictor.bin
```
### Start with the annotation using database constructed above
```bash
vcf_dir="/mnt/Knives/1_data/Sorghum_SNP"
snpeff_dir="/home/liangyu/3_software/snpEff"

### Treat as the entire genome for annotation
cd $snpeff_dir
nohup java -Xmx64g -jar $snpeff_dir/snpEff.jar \
    Sbicolor ${vcf_dir}/filter_mac1_select > ${vcf_dir}/filter_mac1_select.ann.vcf & 

mv $snpeff_dir/snpEff_genes.txt ${vcf_dir}/filter_mac1_snpEff_genes.txt
mv $snpeff_dir/snpEff_summary.html ${vcf_dir}/filter_mac1_snpEff_summary.html
```

## **Pan-genome consturction of the six accessions**
### Merge sub-reads for six accessions
```bash
work_dir="/mnt/Knives/1_data/Sorghum_Pangenome"

### A acession
cat ERR8097262_1.fastq.gz ERR8097502_1.fastq.gz ERR8097508_1.fastq.gz ERR8097613_1.fastq.gz ERR8097641_1.fastq.gz ERR8106149_1.fastq.gz	ERR8106214_1.fastq.gz > $work_dir/1_Reads/PI533871_R1.fastq.gz
cat ERR8097262_2.fastq.gz ERR8097502_2.fastq.gz	ERR8097508_2.fastq.gz ERR8097613_2.fastq.gz ERR8097641_2.fastq.gz ERR8106149_2.fastq.gz	ERR8106214_2.fastq.gz > $work_dir/1_Reads/PI533871_R2.fastq.gz

### B accession
cat ERR8017558_1.fastq.gz ERR8056536_1.fastq.gz ERR8058553_1.fastq.gz ERR8060479_1.fastq.gz ERR8060670_1.fastq.gz ERR8093059_1.fastq.gz ERR8156660_1.fastq.gz ERR8157009_1.fastq.gz ERR8165680_1.fastq.gz ERR8170065_1.fastq.gz > $work_dir/1_Reads/PI533961_R1.fastq.gz
cat ERR8017558_2.fastq.gz ERR8056536_2.fastq.gz ERR8058553_2.fastq.gz ERR8060479_2.fastq.gz ERR8060670_2.fastq.gz ERR8093059_2.fastq.gz	ERR8156660_2.fastq.gz ERR8157009_2.fastq.gz ERR8165680_2.fastq.gz ERR8170065_2.fastq.gz > $work_dir/1_Reads/PI533961_R2.fastq.gz

### C accession
cat ERR8097252_1.fastq.gz ERR8097404_1.fastq.gz ERR8097551_1.fastq.gz ERR8097625_1.fastq.gz ERR8106246_1.fastq.gz ERR8157191_1.fastq.gz ERR8281602_1.fastq.gz > $work_dir/1_Reads/PI656041_R1.fastq.gz
cat ERR8097252_2.fastq.gz ERR8097404_2.fastq.gz ERR8097551_2.fastq.gz ERR8097625_2.fastq.gz ERR8106246_2.fastq.gz ERR8157191_2.fastq.gz ERR8281602_2.fastq.gz > $work_dir/1_Reads/PI656041_R2.fastq.gz

### D accession
cat ERR8016753_1.fastq.gz ERR8056538_1.fastq.gz ERR8058451_1.fastq.gz ERR8060209_1.fastq.gz ERR8093075_1.fastq.gz ERR8097425_1.fastq.gz ERR8156910_1.fastq.gz ERR8169608_1.fastq.gz ERR8169624_1.fastq.gz ERR8281567_1.fastq.gz > $work_dir/1_Reads/PI656053_R1.fastq.gz
cat ERR8016753_2.fastq.gz ERR8056538_2.fastq.gz ERR8058451_2.fastq.gz ERR8060209_2.fastq.gz ERR8093075_2.fastq.gz ERR8097425_2.fastq.gz ERR8156910_2.fastq.gz ERR8169608_2.fastq.gz ERR8169624_2.fastq.gz ERR8281567_2.fastq.gz > $work_dir/1_Reads/PI656053_R2.fastq.gz

### E accession
cat ERR8049945_1.fastq.gz ERR8060163_1.fastq.gz ERR8093118_1.fastq.gz ERR8093131_1.fastq.gz ERR8157017_1.fastq.gz ERR8168867_1.fastq.gz ERR8172437_1.fastq.gz ERR8281442_1.fastq.gz ERR8281537_1.fastq.gz ERR8281553_1.fastq.gz > $work_dir/1_Reads/PI656076_R1.fastq.gz
cat ERR8049945_2.fastq.gz ERR8060163_2.fastq.gz ERR8093118_2.fastq.gz ERR8093131_2.fastq.gz ERR8157017_2.fastq.gz ERR8168867_2.fastq.gz ERR8172437_2.fastq.gz ERR8281442_2.fastq.gz ERR8281537_2.fastq.gz ERR8281553_2.fastq.gz > $work_dir/1_Reads/PI656076_R2.fastq.gz

### F accession
cat ERR8017472_1.fastq.gz ERR8056276_1.fastq.gz ERR8060401_1.fastq.gz ERR8156984_1.fastq.gz ERR8168744_1.fastq.gz ERR8168865_1.fastq.gz ERR8168964_1.fastq.gz ERR8169014_1.fastq.gz ERR8169497_1.fastq.gz ERR8169853_1.fastq.gz > $work_dir/1_Reads/PI656096_R1.fastq.gz
cat ERR8017472_2.fastq.gz ERR8056276_2.fastq.gz ERR8060401_2.fastq.gz ERR8156984_2.fastq.gz ERR8168744_2.fastq.gz ERR8168865_2.fastq.gz ERR8168964_2.fastq.gz ERR8169014_2.fastq.gz ERR8169497_2.fastq.gz ERR8169853_2.fastq.gz > $work_dir/1_Reads/PI656096_R2.fastq.gz

### Check the reads numbers of merged fastaq files
for i in $(ls $work_dir/1_Reads/*.fastq.gz);
do
	echo $(zcat $i|wc -l)/4|bc >> $work_dir/1_Reads/combined_read.counts
done
```

***Check numbers of reads after merging
```
Row Labels	Sum of read_counts
PI533871_1	112181186
PI533871_2	112181186
PI533961_1	88840410
PI533961_2	88840410
PI656041_1	128287213
PI656041_2	128287213
PI656053_1	90037052
PI656053_2	90037052
PI656076_1	88640725
PI656076_2	88640725
PI656096_1	92877627
PI656096_2	92877627
```
### Trim adaptor sequences (Treat as paired end)
```bash
raw_reads_dir="/mnt/Knives/1_data/Sorghum_Pangenome/1_Reads"
trim_reads_dir="/mnt/Knives/1_data/Sorghum_Pangenome/1_Reads/Clean_reads"

for i in PI533871 PI533961 PI656041 PI656053 PI656076 PI656096; 
do
 trimmomatic PE -threads 10 \
    ${raw_reads_dir}/${i}_1.fastq.gz \
    ${raw_reads_dir}/${i}_2.fastq.gz \
    ${trim_reads_dir}/${i}_forward_paired.fq.gz ${trim_reads_dir}/${i}_forward_unpaired.fq.gz \
    ${trim_reads_dir}/${i}_reverse_paired.fq.gz ${trim_reads_dir}/${i}_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/home/liangyu/anaconda3/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2 LEADING:3 TRAILING:3 MINLEN:36

   ### check the numbers of reads within the unpaired reads file
   echo $(zcat ${trim_reads_dir}/${i}_forward_paired.fq.gz | wc -l)/4|bc >> ${trim_reads_dir}/forward_clean.count
   echo $(zcat ${trim_reads_dir}/${i}_reverse_paired.fq.gz | wc -l)/4|bc >> ${trim_reads_dir}/reverse_clean.count

   echo $(zcat ${trim_reads_dir}/${i}_forward_unpaired.fq.gz | wc -l)/4|bc >> ${trim_reads_dir}/forward_unpair.count
   echo $(zcat ${trim_reads_dir}/${i}_reverse_unpaired.fq.gz | wc -l)/4|bc >> ${trim_reads_dir}/reverse_unpair.count

done

### move all untrimmed reads as archives
mv /mnt/Knives/1_data/Sorghum_Pangenome/1_Reads/*.fastq.gz /mnt/Kanta/Sorghum_reads/
```
### Genome asembly using Megahit pipeline
The manual of megahit can be refered [**here**](https://github.com/voutcn/megahit)
```bash
trim_reads_dir="/mnt/Knives/1_data/Sorghum_Pangenome/1_Reads/Clean_reads"
assembly_dir="/mnt/Knives/1_data/Sorghum_Pangenome/2_Assembly"

# Genome assembly using megahit
for i in PI533871 PI533961 PI656041 PI656053 PI656076 PI656096; 
do
    mkdir ${assembly_dir}/${i}
    megahit \
        -t 10 \
        -1 ${trim_reads_dir}/${i}_forward_paired.fq.gz \
        -2 ${trim_reads_dir}/${i}_reverse_paired.fq.gz \
        -o ${assembly_dir}/${i}/$i \
        2>$i/$i.log

    # Remove unnecessary files
    rm ${assembly_dir}/${i}/$i/checkpoints.txt \
        ${assembly_dir}/${i}/$i/done \
        ${assembly_dir}/${i}/$i/log \
        ${assembly_dir}/${i}/$i/options.json

    rm -r ${assembly_dir}/${i}/$i/intermediate_contigs
done
```
### Evaluation of genome assembly by BUSCO
[**BUSCO**](https://busco.ezlab.org/busco_userguide.html#conda-package) maunal can be referred 

```bash
BUSCO_dir="/mnt/Knives/1_data/Sorghum_Pangenome/6_BUSCO"
assembly_dir="/mnt/Knives/1_data/Sorghum_Pangenome/2_Assembly"

### Check the databse
busco --list-datasets

### Genome assembly using megahit
for i in PI533871 PI533961 PI656041 PI656053 PI656076 PI656096; 
do
    ln -s  ${assembly_dir}/${i}/${i}/final.contigs.fa ${BUSCO_dir}/${i}.fasta

    ### start the estimation
    busco \
    -i ${BUSCO_dir}/${i}.fasta \
    -l embryophyta_odb10 \
    -o ${i} \
    -m genome \
    --cpu 8 \
    -f
done
```


### Evaluation of genome assembly by QUAST
[**QUAST**](http://quast.sourceforge.net/docs/manual.html#sec1) Manual can be refered
```bash
# Working directory
quast_dir="/mnt/Knives/1_data/Sorghum_Pangenome/3_Quast"
megahit_dir="/mnt/Knives/1_data/Sorghum_Pangenome/2_megahit"

for i in $(cat $quast_dir/Sbicolor.list);
do 
    quast \
        -o ${quast_dir}/${i}_quast_identity90 \
        ${megahit_dir}/${i}/final.contigs.fa \
        -r ${quast_dir}/Sbicolor.fa \
        --threads 10 \
        --min-identity 90 \
        --min-alignment 500 \
        --fragmented \
        --no-html \
        --no-icarus \
        --no-snps \
        --eukaryote
done 
```
**Summary of the assembled sequences**
```bash
# Assembled genome size
for i in $(cat $quast_dir/Sbicolor.list); 
do 
    cat {quast_dir}/${i}_quast_identity90/report.txt | sed -n '10p' | awk '{print $6}' >> {quast_dir}/ASSEMBLED.size
done

# Fully unaligned sequence length
for i in $(cat $quast_dir/Sbicolor.list); 
do 
    cat {quast_dir}/${i}_quast_identity90/contigs_reports/unaligned_report.txt | sed -n '5p' | awk '{print $4}' >> {quast_dir}/Unaligned.size
done

# Partially unaligned length
for i in $(cat $quast_dir/Sbicolor.list);
do 
    cat {quast_dir}/${i}_quast_identity90/contigs_reports/unaligned_report.txt | sed -n '7p' | awk '{print $4}' >> {quast_dir}/P_aligned.size
done
```

### Extract unaligned sequences
[**Scripts**](https://github.com/Sunhh/NGS_data_processing) for cleaning sequences can be accessed under the Github page of Sunhh
List of perl scripts used in studies
```classify_region_byBn6.pl
fasta-len.pl
fasta_change_line.pl
get_Ex_region.pl
get_ctg.pl
parse_bn6.pl
remove_seq_from_id.pl
```
```bash
scripts_dir="/mnt/Knives/1_data/Sorghum_Pangenome/3_Quast/1_scripts"

perl $scripts_dir/get_ctg.pl $quast_dir/Sbicolor.list
# In this perl script, it will read the accessions in the given list and extract the unaligned sequences.
### Parameters (sequences with less than 400bp was dropped)
### remove the inter-mediate files
rm *changedLine.fasta *changedLine.fasta.index.dir *changedLine.fasta.index.pag

### Check the integrated unaligned sequences
grep -v ">" PI533871.unalign.fasta | tr -d '\n' | wc -c
grep -v ">" PI533961.unalign.fasta | tr -d '\n' | wc -c
grep -v ">" PI656041.unalign.fasta | tr -d '\n' | wc -c
grep -v ">" PI656053.unalign.fasta | tr -d '\n' | wc -c
grep -v ">" PI656076.unalign.fasta | tr -d '\n' | wc -c
grep -v ">" PI656096.unalign.fasta | tr -d '\n' | wc -c

gzip *unalign.fasta

# Combine the unaligned sequences
zcat $quast_dir/*unalign.fasta.gz > $quast_dir/S.boicolor.unalign.fa
```

### Remove duplicates
**PLEASE refer manual and scripts uner the Honghe Sun's [Github](https://github.com/Sunhh/NGS_data_processing) page**

The purpose of this step is to remove redundant and duplicated sequences from unaligned sequences derived from differnet accessios. 
```bash
### set directory
blast_dir="/mnt/Knives/1_data/Sorghum_Pangenome/3_Quast/2_Blast"
scripts_dir="/mnt/Knives/1_data/Sorghum_Pangenome/3_Quast/1_scripts"
quast_dir="/mnt/Knives/1_data/Sorghum_Pangenome/3_Quast"

### Merge the unaligned sequences and reference genome as blastndb 
cat $quast_dir/S.boicolor.unalign.fa $quast_dir/Sbicolor.fa > $quast_dir/Sbicolor_pan.fa

### Build the blastdb
makeblastdb -in $quast_dir/Sbicolor_pan.fa \
    -dbtype nucl \
    -out $blast_dir/Sbicolor_Panref \
    -parse_seqids

### Run blastn with identity of 90 as cut-off
blastn -query $quast_dir/S.boicolor.unalign.fa \
    -out $blast_dir/Sbicolor_unalign.out\
    -db  $blast_dir/Sbicolor_Panref \
    -word_size 25 -perc_identity 90 -dust no -num_threads 20 -max_target_seqs 10 -evalue 1e-10 \
    -outfmt '7 qseqid sseqid qlen slen length qstart qend sstart send sstrand pident nident mismatch gapopen gaps qcovhsp qcovs score bitscore evalue' 


# Remove redundancy from blast
python2 $scripts_dir/Remove_redundancy_by_blast.py \
    $quast_dir/S.boicolor.unalign.fa \
    $blast_dir/Sbicolor_unalign.out 90 300 \
    > $blast_dir/S.boicolor.unalign.clean.fa

# Check the number of sequences
grep ">" $blast_dir/S.boicolor.unalign.clean.fa | wc

# Total length of sequences
grep -v ">" $blast_dir/S.boicolor.unalign.clean.fa | tr -d '\n' | wc -c 
```
### Remove contamination
#### Download the NCBI nr/nt database
The link for [**NCBI databse**](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/)
```bash
Database_dir="/mnt/Knives/1_data/Sorghum_Pangenome/0_data"

### Start with downloading database from NT database
ascp -k 1 -QT -l 500m \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/FASTA/nt.gz \
    $Database_dir/

gunzip $Database_dir/nt.gz

### create database index
makeblastdb -in $Database_dir/nt \
    -dbtype nucl \
    -out $Database_dir/nt \
    -parse_seqids

### Download the formated databse in the creation failed 
for i in $(seq -w 00 70);
do
	ascp -k 1 -QT -l 500m \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/nt.$i.tar.gz \
    $Database_dir/

    tar -zxvf nt.$i.tar.gz
done

```
[**Trouble shotting**](https://www.cxymm.net/article/zhangjunya/108279859) information of blastndb version was attached 
The blast+ should be unudated. Download and intall the most recent [**verison**](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 

```bash
 wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz

 ### add the new blast to teh path 
echo 'export PATH=/home/liangyu/ncbi-blast-2.13.0+/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```


#### Perform the cleaning steps
```bash
# Working directory
Con_dir="/mnt/Knives/1_data/Sorghum_Pangenome/4_RemoveContamination"
quast_dir="/mnt/Knives/1_data/Sorghum_Pangenome/3_Quast"
scripts_dir="/mnt/Knives/1_data/Sorghum_Pangenome/3_Quast/1_scripts"
blast_dir="/mnt/Knives/1_data/Sorghum_Pangenome/3_Quast/2_Blast"

# Run blastn using the NCBI database as the subject sequences
blastn -task megablast \
    -num_threads 24 \
    -max_target_seqs 20 \
    -evalue 1e-5 \
    -dust yes \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sstrand qsscinames sskingdoms stitle" \
    -db /mnt/Knives/1_data/Sorghum_Pangenome/0_data/nt \
    -query ${blast_dir}/S.boicolor.unalign.clean.fa \
    -out $Con_dir/S.boicolor.unalign.bn6

# parse blastn results
perl $scripts_dir/classify_region_byBn6.pl \
    $Con_dir/S.boicolor.unalign.bn6 \
    -maxUn 50 \
    -InList 'Eukaryota:Satellite:rDNA' \
    -ExList 'NA:Viruses:Bacteria:Archaea:Chloroplast:Mitochondrion:Plastid' \
    -joinInEx $Con_dir/S.boicolor.unalign.blastn.ExJn \
    > $Con_dir/S.boicolor.unalign.ExJn.list \
    2> $Con_dir/S.boicolor.unalign.ExJn.log

perl $scripts_dir/get_Ex_region.pl \
    $Con_dir/S.boicolor.unalign.blastn.ExJn \
    1> $Con_dir/S.boicolor.unalign.ExSeq 2> $Con_dir/S.boicolor.unalign.ExID

# Extract unaligned clean sequences
perl  $scripts_dir/remove_seq_from_id.pl \
    $Con_dir/S.boicolor.unalign.ExID \
    $blast_dir/S.boicolor.unalign.clean.fa \
    $Con_dir/S.boicolor.RemoveContanmination.fa

# Check the number of sequences
grep ">" $Con_dir/S.boicolor.RemoveContanmination.fa | wc
# Total length of sequences
grep -v ">" $Con_dir/S.boicolor.RemoveContanmination.fa | tr -d '\n' | wc -c 
```

### Novel gene annotation and prediciton
#### Using lift-off mapping strategy
[**lift-off**](https://github.com/agshumate/Liftoff) protol can be read from github page 

```bash
```

#### Using RNA-seq reads mapping strategy
 
#### Mapping reads to Clean Pan sequences
```bash
fastq_dir="/home/liangyu/1_Project/8_Sorghum/1_Salmon/0_data"
Ref_dir="/mnt/Knives/1_data/Sorghum_Pangenome/4_RemoveContamination"
sample_dir="/home/liangyu/1_Project/8_Sorghum/1_Salmon/0_data"
output_dir="/home/liangyu/1_Project/8_Sorghum/1_Salmon-Pan"

###build index for salmon
/home/liangyu/anaconda3/envs/salmon/bin/salmon index -t $Ref_dir/S.boicolor.RemoveContanmination.fa -i $Ref_dir/S.boicolor.RemoveContanmination.salmon

IFS=$'\n';
for LINE in $(cat ${sample_dir}/sample.list);
do

	fq_file=$(echo ${LINE} | awk '{ print $1}')
	salmon_file=$(echo ${LINE} | awk '{ print $2}')

###perform quantification for each sample
	/home/liangyu/anaconda3/envs/salmon/bin/salmon quant -l A -1 $fastq_dir/${fq_file}_R1.fastq.gz -2 $fastq_dir/${fq_file}_R2.fastq.gz -i $Ref_dir/S.boicolor.RemoveContanmination.salmon -o ${output_dir}/${salmon_file}.salmon.count -p 10
done
```


# **Section II Gene expression and regulatory network analysis of experiments**
## STEP 1 Quantification of reads counts using featureCounts 
The step is to count the reads per primary transcripts using [**featureCounts**](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html). 

### **Extraction of the primary transcripts and longest transcripts for lncRNAs based on sorghum gene annotation**
```bash
### for coding genes
gff_dir="/home/liangyu/1_Project/3_RNA-mods/8_gffcompare"
### Generate a total of 34129 out of 47121  primary lncRNAs as representative 
grep '>' $gff_dir/Sorghum_cds_primary_TrasncriptsOnly.fasta > $gff_dir/Sb_PT.list 
for i in $(cat $gff_dir/Sb_PT.list);
do
	grep -w $i $gff_dir/Sbicolor_454_v3.1.1.gene.gff3 >> $gff_dir/Sorghum.gff_primary_transcripts.gff3
done

```
**Example header of Sorghum.gff_primary_transcripts.gff3**\
exon as "feature"\
Parent as "meta-fature"/"attribute"
```
Chr10	phytozomev12	mRNA	6867058	6894701	.	+	.	ID=Sobic.010G080700.2;Name=Sobic.010G080700.2	pacid=37908434	longest=1	ancestorIdentifier=Sobic.010G080700.2.v2.1
Chr10	phytozomev12	five_prime_UTR.	6867058	6867546	.	+	.	ID=Sobic.010G080700.2.five_prime_UTR.1;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	five_prime_UTR.	6867920	6868028	.	+	.	ID=Sobic.010G080700.2.five_prime_UTR.2;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6868029	6868060	.	+	0	ID=Sobic.010G080700.2:exon:1;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6868136	6868238	.	+	1	ID=Sobic.010G080700.2:exon:2;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6868579	6868700	.	+	0	ID=Sobic.010G080700.2:exon:3;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6868787	6868879	.	+	1	ID=Sobic.010G080700.2:exon:4;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6870395	6870461	.	+	1	ID=Sobic.010G080700.2:exon:5;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6870544	6870634	.	+	0	ID=Sobic.010G080700.2:exon:6;Parent=Sobic.010G080700.2	pacid=37908434		
Chr10	phytozomev12	exon	6870733	6870785	.	+	2	ID=Sobic.010G080700.2:exon:7;Parent=Sobic.010G080700.2	pacid=37908434		

```
```bash
### for lncRNAs
gff_dir="/home/liangyu/1_Project/3_RNA-mods/8_gffcompare"
### Generate a total of 79054 out of 84720  primary lncRNAs as representative 
for i in $(cat $gff_dir/Sorghum_lncRNA.list);
do
	grep -w $i $gff_dir/Sorghum_lncRNAs.gtf >> $gff_dir/Sorghum_primary_lncRNAs.gtf
done

### Check the numbers of lncRNAs in primary_lncRNAs gtf file (79054)
grep "transcript" $gff_dir/Sorghum_primary_lncRNAs.gtf | wc -l 
```
**Example header of Sorghum_primary_lncRNAs.gtf**\
exon as "feature"\
geneID as "meta-fature"/"attribute"
```
Chr1    Evolinc transcript      266228  267820  .       +       .       ID=TCONS_00000001;geneID=XLOC_000001
Chr1    Evolinc exon    266228  267820  .       +       .       ID=TCONS_00000001
Chr1    Evolinc transcript      330591  334048  .       +       .       ID=TCONS_00000002;geneID=XLOC_000002
Chr1    Evolinc exon    330591  330746  .       +       .       ID=TCONS_00000002
Chr1    Evolinc exon    333896  334048  .       +       .       ID=TCONS_00000002
Chr1    Evolinc transcript      786661  787319  .       +       .       ID=TCONS_00000003;geneID=XLOC_000003
Chr1    Evolinc exon    786661  786977  .       +       .       ID=TCONS_00000003
Chr1    Evolinc exon    787147  787319  .       +       .       ID=TCONS_00000003
Chr1    Evolinc transcript      947283  947986  .       +       .       ID=TCONS_00000004;geneID=XLOC_000004
Chr1    Evolinc exon    947283  947986  .       +       .       ID=TCONS_00000004
Chr1    Evolinc transcript      1018713 1020325 .       +       .       ID=TCONS_00000005;geneID=XLOC_000005
Chr1    Evolinc exon    1018713 1019320 .       +       .       ID=TCONS_00000005
Chr1    Evolinc exon    1019387 1019751 .       +       .       ID=TCONS_00000005
Chr1    Evolinc exon    1019900 1020325 .       +       .       ID=TCONS_00000005

```

### **Perform the featureCounts for coding genes at exon level**
```bash
gff_dir="/home/liangyu/1_Project/3_RNA-mods/8_gffcompare"
BAM_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/1_Sample" 
output_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/6_FeatureCounts_results"

###Perform featureCounts for coding genes
featureCounts -T 20 -s 1 -p \
  -a $GTF_dir/Sorghum.gff_primary_transcripts.gff3 \
  -o $output_dir/Sorghum_transcript.count \
  -t exon \
  -g Parent\
  -O \
  $BAM_dir/*.bam > transcripts_count.log
  
grep "Total reads" $output_dir/transcript_exon_count.log | cut -d " " -f8 > $output_dir/transcript_exon_Total_reads

###Clean GeneID
grep -v 'featureCounts v1.6.0; Command' $output_dir/Sorghum_transcript.count | cut -d ":" -f2 > Sorghum_transcript_clean.count
sed -i 's@/mnt/Knives/1_data/5_MODs_bam/sorghum/@@g' $output_dir/Sorghum_transcript_clean.count
sed -i "s@/mnt/Knives/1_data/5_MODs_bam/sorghum/1_Sample/@@g" Sorghum_transcript.count.summary
```
### check if meta-feature and features counted correctly for transcripts
```
Features : 154048                                                       
Meta-features : 34129                                                   
Chromosomes/contigs : 90
```
### **Perform the featureCounts for lncRNAs (exon level)**
```bash
gff_dir="/home/liangyu/1_Project/3_RNA-mods/8_gffcompare"
BAM_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/1_Sample" 
output_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/6_FeatureCounts_results"

###Perform featureCounts for lncRNAs
nohup featureCounts -T 20 -s 1 -p \
  -a $GTF_dir/Sorghum.gff_primary_transcripts.gff3 \
  -o $output_dir/Sorghum_lncRNA.count \
  -t exon \
  -g transcript \
  -O \
  $BAM_dir/*.bam > lncRNA_exon_count.log 

###Clean GeneID
grep -v 'featureCounts v1.6.0; Command' $output_dir/Sorghum_lncRNA.count | cut -d ":" -f2 > Sorghum_lncRNA_clean.count 
sed -i 's@/mnt/Knives/1_data/5_MODs_bam/sorghum/@@g' $output_dir/Sorghum_lncRNA_clean.count
sed -i "s@/mnt/Knives/1_data/5_MODs_bam/sorghum/1_Sample/@@g" Sorghum_lncRNA.count.summary
```
### check if meta-feature and features counted correctly for lncRNAs
```
Features : 87989                                                    
Meta-features : 79054   (grouped features)                                                
Chromosomes/contigs : 10  
```

### **Perform the featureCounts for modififed sites (per PTM sites level)**
```bash
gff_dir="/home/liangyu/1_Project/3_RNA-mods/8_gffcompare"
BAM_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/1_Sample" 
output_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/6_FeatureCounts_results"

###Perform featureCounts for coding genes
featureCounts -T 20 -s 1 -p \
  -a $gff_dir/Sorghum_PTM.gff3 \
  -o $output_dir/Sorghum_PTM.count \
  -t exon \
  -g ID\
  -O \
  $BAM_dir/*.bam > PTM_count.log
  
grep "Total reads" $output_dir/PTM_count.log | cut -d " " -f8 > $output_dir/PTM_Total_reads

###Clean GeneID
grep -v 'featureCounts v1.6.0; Command' $output_dir/Sorghum_PTM.count | cut -d ":" -f2 > Sorghum_PTM_clean.count
sed -i 's@/mnt/Knives/1_data/5_MODs_bam/sorghum/@@g' $output_dir/Sorghum_PTM_clean.count
sed -i "s@/mnt/Knives/1_data/5_MODs_bam/sorghum/1_Sample/@@g" Sorghum_PTM.count.summary
```

### check if meta-feature and features counted correctly for PTMs
```
Features : 20126                                                        
Meta-features : 20126                                                   
Chromosomes/contigs : 11 
```

## STEP 2 Normalization of read counts and differentially expressed genes (DEGs) analysis
### **Install or load packages** 
```r
library(DESeq2)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(knitr)
library(RColorBrewer)
library(corrplot)
library(PoiClaClu)
library(plotly)
library(stats)

###NOTE: This step will be excuted only if you've saved working env. Or should just skip
load("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/2_Pipelines/DEseq2_FeatureCounts.RData")
```
### **Import FeatureCounts data**
```r
#Input list of genes from one dataset been modififed 
FC1 <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/3_FeatureCounts/Sorghum_20220310_FC_final.csv", header = T)
info <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/3_FeatureCounts/Sorghum_20220310-info_final.csv", header = T)

#Place Gene_IDs into appropriate format for DESeq2
row.names(FC1) <- FC1[,1]
#FC <- FC1[,3:102]
FC <- FC1[,3:98]

info$Treatment <- as.factor(info$Treatment)
rownames(info) <- info$Sample

info$Group <- paste0(info$Accession, "_", info$Treatment, "_", info$Time)

### Test of all rownames from info (Coldata) matched the colnames in FC (counts_data)
all(colnames(FC) %in% rownames(info))

### Test of all rownames from info (Coldata) matched the order of colnames in FC (counts_data)
all(colnames(FC) == rownames(info))
```
### **Experimental design and DEGs analysis**
```r
#start normalization
dds <- DESeqDataSetFromMatrix(countData = FC,
                              colData = info, 
                              design = ~ Group)

## Keep the lines with great than 1 count (at least one count per gene)
keep <- rownames(dds)[rowSums(counts(dds)) > 96]
dds <- dds[keep,]

### check numbers of genes remained for analysis
nrow(dds)
###Relevel the factor
dds$Treatment <- relevel(dds$Treatment, ref = "WW")

### Reads normalization and DEGs analysis
dds <- DESeq(dds)

### Export data for A accessions
resA1 <- as.data.frame(results(dds, contrast=c("Group","A_WL_TP1","A_WW_TP1")))
resA1$Sample = "A_TP1"
resA3 <- as.data.frame(results(dds, contrast=c("Group","A_WL_TP3","A_WW_TP3")))
resA3$Sample = "A_TP3"
resA5 <- as.data.frame(results(dds, contrast=c("Group","A_WL_TP5","A_WW_TP5")))
resA5$Sample = "A_TP5"
resA7 <- as.data.frame(results(dds, contrast=c("Group","A_WL_TP7","A_WW_TP7")))
resA7$Sample = "A_TP7"
write.csv(rbind(resA1, resA3, resA5, resA7), "C:/Users/Leon/Desktop/DEGs_A.csv")

### Export data for B accessions
resB1 <- as.data.frame(results(dds, contrast=c("Group","B_WL_TP1","B_WW_TP1")))
resB1$Sample = "B_TP1"
resB3 <- as.data.frame(results(dds, contrast=c("Group","B_WL_TP3","B_WW_TP3")))
resB3$Sample = "B_TP3"
resB5 <- as.data.frame(results(dds, contrast=c("Group","B_WL_TP5","B_WW_TP5")))
resB5$Sample = "B_TP5"
resB7 <- as.data.frame(results(dds, contrast=c("Group","B_WL_TP7","B_WW_TP7")))
resB7$Sample = "B_TP7"
write.csv(rbind(resB1, resB3, resB5, resB7), "C:/Users/Leon/Desktop/DEGs_B.csv")

### Export data for C accessions
resC1 <- as.data.frame(results(dds, contrast=c("Group","C_WL_TP1","C_WW_TP1")))
resC1$Sample = "C_TP1"
resC3 <- as.data.frame(results(dds, contrast=c("Group","C_WL_TP3","C_WW_TP3")))
resC3$Sample = "C_TP3"
resC5 <- as.data.frame(results(dds, contrast=c("Group","C_WL_TP5","C_WW_TP5")))
resC5$Sample = "C_TP5"
resC7 <- as.data.frame(results(dds, contrast=c("Group","C_WL_TP7","C_WW_TP7")))
resC7$Sample = "C_TP7"
write.csv(rbind(resC1, resC3, resC5, resC7), "C:/Users/Leon/Desktop/DEGs_C.csv")

### Export data for D accessions
resD1 <- as.data.frame(results(dds, contrast=c("Group","D_WL_TP1","D_WW_TP1")))
resD1$Sample = "D_TP1"
resD3 <- as.data.frame(results(dds, contrast=c("Group","D_WL_TP3","D_WW_TP3")))
resD3$Sample = "D_TP3"
resD5 <- as.data.frame(results(dds, contrast=c("Group","D_WL_TP5","D_WW_TP5")))
resD5$Sample = "D_TP5"
resD7 <- as.data.frame(results(dds, contrast=c("Group","D_WL_TP7","D_WW_TP7")))
resD7$Sample = "D_TP7"
write.csv(rbind(resD1, resD3, resD5, resD7), "C:/Users/Leon/Desktop/DEGs_D.csv")

### Export data for E accessions
resE1 <- as.data.frame(results(dds, contrast=c("Group","E_WL_TP1","E_WW_TP1")))
resE1$Sample = "E_TP1"
resE3 <- as.data.frame(results(dds, contrast=c("Group","E_WL_TP3","E_WW_TP3")))
resE3$Sample = "E_TP3"
resE5 <- as.data.frame(results(dds, contrast=c("Group","E_WL_TP5","E_WW_TP5")))
resE5$Sample = "E_TP5"
resE7 <- as.data.frame(results(dds, contrast=c("Group","E_WL_TP7","E_WW_TP7")))
resE7$Sample = "E_TP7"
write.csv(rbind(resE1, resE3, resE5, resE7), "C:/Users/Leon/Desktop/DEGs_E.csv")

### Export data for F accessions
resF1 <- as.data.frame(results(dds, contrast=c("Group","F_WL_TP1","F_WW_TP1")))
resF1$Sample = "F_TP1"
resF3 <- as.data.frame(results(dds, contrast=c("Group","F_WL_TP3","F_WW_TP3")))
resF3$Sample = "F_TP3"
resF5 <- as.data.frame(results(dds, contrast=c("Group","F_WL_TP5","F_WW_TP5")))
resF5$Sample = "F_TP5"
resF7 <- as.data.frame(results(dds, contrast=c("Group","F_WL_TP7","F_WW_TP7")))
resF7$Sample = "F_TP7"
write.csv(rbind(resF1, resF3, resF5, resF7), "C:/Users/Leon/Desktop/DEGs_F.csv")
```
### **Export data for inter-accession comparison**
**(Using the C accession as reference accession for both conditions)**
For example, for a simple two-group comparison, this would return the log2 fold changes of the second group over the first group (the reference level). Please see examples below and in the vignette.

```r
PATH = "C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/10_DEGs/Inter-accessions"
resAC_WL <-  as.data.frame(results(dds, contrast=c("Group","C_WL_TP7","A_WL_TP7")))
write.csv(resAC_WL, paste0(PATH, "/DEGs_resAC_WL.csv"))
resAC_WW <-  as.data.frame(results(dds, contrast=c("Group","C_WW_TP7","A_WW_TP7")))
write.csv(resAC_WW, paste0(PATH, "/DEGs_resAC_WW.csv"))

resBC_WL <-  as.data.frame(results(dds, contrast=c("Group","C_WL_TP7","B_WL_TP7")))
write.csv(resBC_WL, paste0(PATH, "/DEGs_resBC_WL.csv"))
resBC_WW <-  as.data.frame(results(dds, contrast=c("Group","C_WW_TP7","B_WW_TP7")))
write.csv(resBC_WW, paste0(PATH, "/DEGs_resBC_WW.csv"))

resDC_WL <-  as.data.frame(results(dds, contrast=c("Group","C_WL_TP7","D_WL_TP7")))
write.csv(resDC_WL, paste0(PATH, "/DEGs_resDC_WL.csv"))
resDC_WW <-  as.data.frame(results(dds, contrast=c("Group","C_WW_TP7","D_WW_TP7")))
write.csv(resDC_WW, paste0(PATH, "/DEGs_resDC_WW.csv"))

resEC_WL <-  as.data.frame(results(dds, contrast=c("Group","C_WL_TP7","E_WL_TP7")))
write.csv(resEC_WL, paste0(PATH, "/DEGs_resEC_WL.csv"))
resEC_WW <-  as.data.frame(results(dds, contrast=c("Group","C_WW_TP7","E_WW_TP7")))
write.csv(resEC_WW, paste0(PATH, "/DEGs_resEC_WW.csv"))

resFC_WL <-  as.data.frame(results(dds, contrast=c("Group","C_WL_TP7","F_WL_TP7")))
write.csv(resFC_WL, paste0(PATH, "/DEGs_resFC_WL.csv"))
resFC_WW <-  as.data.frame(results(dds, contrast=c("Group","C_WW_TP7","F_WW_TP7")))
write.csv(resFC_WW, paste0(PATH, "/DEGs_resFC_WW.csv"))
```
### **Extraction of normalized counts**
```r
#Export the count information
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- data.frame(normalized_counts)

normalized_counts$Geneid <- rownames(normalized_counts)
normalized_counts <- left_join(normalized_counts, FC1[,c(1,2)], by = "Geneid")
#NFC  <- normalized_counts[,c(97,98,1:96)]
NFC  <- normalized_counts[,c(101,102,1:100)]

### 30248 genes retained for analysis
write.table(NFC, file="C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/3_FeatureCounts/Sorghum_update.csv", sep="\t", quote=F, col.names=NA)
```
### **PCA based on vst transfomration**
```r
### Use vsd transformation
vsd <- vst(dds, blind = F)

### Modified version of PCA to extract PC1, PC2, and PC3
plotPCA <- function (object, intgroup=c("Accession", "Treatment"), ntop = 2332, returnData=TRUE) 
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:4]
        return(d)
    }
    ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC3: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC4: ", round(percentVar[2] * 
        100), "% variance")) + coord_fixed()
}

### extract the PC1, PC2, PC3, and PC4 data
pcaData <- plotPCA(vsd, intgroup=c("Accession", "Time"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar

### Plot PC1 and PC2 by Accessions and treatment
ggplot(pcaData, aes(PC1, PC2, color=Accession, shape=Time)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(2,16))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw()


### Export PCA values for each component
write.csv(pcaData, "C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/5_Results/vsd_PCA_top10SD_100.csv")

#Subset different Accessions for independent WGCNA 
vsd_sub1 <- vsd[ , vsd$Accession == "A" ]
plotPCA(vsd_sub1, intgroup=c("Time"))

vsd_sub2 <- vsd[ , vsd$Accession == "B" ]
plotPCA(vsd_sub2, intgroup=c("Time"))

vsd_sub3 <- vsd[ , vsd$Accession == "C" ]
plotPCA(vsd_sub3, intgroup=c("Time"))

vsd_sub4 <- vsd[ , vsd$Accession == "D" ]
plotPCA(vsd_sub4, intgroup=c("Time"))

vsd_sub5 <- vsd[ , vsd$Accession == "E" ]
plotPCA(vsd_sub5, intgroup=c("Time"))

vsd_sub6 <- vsd[ , vsd$Accession == "F" ]
plotPCA(vsd_sub6, intgroup=c("Time"))
```
### **Distance Matrix clustering under WL and WW condition**
```r
vsd_WL <- vsd[, vsd$Treatment == "WL"]
vsd_WW <- vsd[, vsd$Treatment == "WW"]

sampleDists <- dist(t(assay(vsd_WL)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_WL$Accession, vsd_WL$Time, sep="_")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(brewer.pal(15, "RdYlBu"))(255)
pheatmap(sampleDistMatrix,color = colors,  
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)

sampleDists <- dist(t(assay(vsd_WW)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_WW$Accession, vsd_WW$Time, sep="_")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(brewer.pal(15, "RdYlBu"))(255)
pheatmap(sampleDistMatrix,color = colors,  
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)
```
### **Transform into TPM**
```r
##Calculate the TPM based on definition 
kb <- NFC$Length / 1000
countdata <- NFC[,3:102]
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)

TPM <- data.frame("Name" = NFC$Geneid, tpm)
write.csv(TPM, "C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/3_FeatureCounts/Sorghum_TPM_20220321.csv", row.names = F, quote=FALSE)
```
### **Pearson correaltion for each two samples based on normalized featureCounts/TPM**
```r
library(Hmisc)
vsd_counts <- data.frame(assay(vsd, normalized=TRUE))
nrow(vsd_counts)

#### Plot correlation
P_sample <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/3_Salmon/Exp_correlation.list.csv", header = T)

# Pearson correlation based on Vsd
for (i in 1:nrow(P_sample)){
  Correlation.df <- vsd_counts[,grepl(P_sample[i,1],colnames(vsd_counts))]
  Correlation.plot <- corrplot.mixed(cor(Correlation.df), order = 'AOE')
                      
  assign(P_sample[i,1], Correlation.df)
  assign(paste0(P_sample[i,1],"_Plot"), Correlation.plot)
}

### Calcualte the P-value effects, score of all correlations
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

df1 <- data.frame(matrix(ncol = 4 , nrow = 48))
colnames(df1) <- c("row", "column", "R", "P-value")

for (i in 1:(nrow(df1))){
  Correlation.df <- vsd_counts[,grepl(P_sample[i,1],colnames(vsd_counts))]
  res <- rcorr(Correlation.df[,1], Correlation.df[,2, ],  type=c("pearson"))
  df1[i,1:4] <- flattenCorrMatrix(res$r, res$P)
}

write.csv(P_sample, "/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/3_FeatureCounts/Pearson-correlation_vsd-20220311.csv", row.names = F)

write.csv(vsd_counts, "/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/3_FeatureCounts/Sorghum_20220310_clean_vsd_update100.csv", quote = F)
```

## STEP 4 WGCNA and regulatory network analysis
### **Optimization of the nework construction**

#### Prepare the PARAMETER_TABLE
We would run multiple-times to generate the optimized parameters for network constriction\
Here's the example of meta-table (PARAMETER_TABLE) used for test 
```
#EXPRESSION_FILE #SELECT_NUMBERS #RATIO
Sorghum_20220629_clean_vsd_update96.csv 18297	vsd82
Sorghum_20220629_clean_vsd_update96.csv	19190	vsd86
Sorghum_20220629_clean_vsd_update96.csv	19636	vsd88
Sorghum_20220629_clean_vsd_update96.csv	20083	vsd90
Sorghum_20220629_clean_vsd_update96.csv	20529	vsd92
Sorghum_20220629_clean_vsd_update96.csv	20975	vsd94
Sorghum_20220629_clean_vsd_update96.csv	21421	vsd96
Sorghum_20220629_clean_vsd_update96.csv	21868	vsd98
Sorghum_20220629_clean_vsd_update96.csv	22310	vsd100
```
#### Prepare the R-code
Here's the R code (WGCNA_template.R) used for networks optimizations
```r
library(WGCNA)
library(cluster)
library(dplyr)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(viridis)
library(matrixStats)
library(ggplot2)
library(ggpubr)

enableWGCNAThreads(nThreads = 30)

### load expression file
TPM <- read.csv("/home/liangyu/1_Project/8_Sorghum/3_WGCNA/EXPRESSION_FILE", header = T, row.names = 1)
#Reorder samples orders by colnames
TPM <- TPM[ , order(names(TPM))]

### Calculate the average by loop
for (i in 1:length(vector)){
  datTrait_ave[grepl(vector[i],rownames(datTrait_ave)),] <- colMeans(datTrait[grepl(vector[i],rownames(datTrait)),])
}

TPM_clean <- data.frame(matrix(ncol = 48, nrow = nrow(TPM)))
rownames(TPM_clean) <- rownames(TPM)
colnames(TPM_clean) <- vector

Sample_list <- unique(str_replace(colnames(TPM_clean),"_WL",""))
#Loop calculation for each accesion average
for (i in 1:length(vector)){
  TPM_clean[,grepl(vector[i],colnames(TPM_clean))] <- rowMeans(TPM[,grepl(vector[i],colnames(TPM))])
}

### STEP 1. Load genes with high MAD score
datExpr <- as.data.frame(t(TPM_clean[order(apply(TPM_clean,1,mad), decreasing = T)[1:SELECT_NUMBERS],]))

### STEP 2. Calculate the soft thresholding power (beta)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#Plot the results:
pdf(file = "/home/liangyu/1_Project/8_Sorghum/3_WGCNA/RATIO/Power-estimation_RATIO.pdf", wi = 12, he = 8)

par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
  
pdf(file = "/home/liangyu/1_Project/8_Sorghum/3_WGCNA/RATIO/Topology_RATIO.pdf", wi = 12, he = 8)
k <- softConnectivity(datE=datExpr,power=sft$powerEstimate) 

par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()

sft$powerEstimate

### STEP3.1. Step-by-step network construction and module detection
enableWGCNAThreads(nThreads = 30)

#Co-expression similarity and adjacency
softPower = sft$powerEstimate;
adjacency = adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);

#set the minimum module size
minModuleSize = 100;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 4, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)

### STEP 3.2 Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

### Export the count data
Module_count <- as.data.frame(table(mergedColors))
write.csv(Module_count, "/home/liangyu/1_Project/8_Sorghum/3_WGCNA/RATIO/MODULE_COUNT_RATIO.csv")
save.image("/mnt/Knives/1_data/Sorghum_WGCNA/Sorghum_RATIOscore.RData")
```
#### Prepare the bash file for loop Rcode
```bash
work_dir="/home/liangyu/1_Project/8_Sorghum/3_WGCNA"
env_dir="/mnt/Knives/1_data/Sorghum_WGCNA"

IFS=$'\n';
for LINE in $(cat $work_dir/PARAMETER_TABLE | sed '1d');do

    EXPRESSION_FILE=$(echo ${LINE} | awk '{ print $1 }')
    SELECT_NUMBERS=$(echo ${LINE} | awk '{ print $2 }')
    RATIO=$(echo ${LINE} | awk '{ print $3 }')

    sed "s/EXPRESSION_FILE/$EXPRESSION_FILE/g" $work_dir/WGCNA_template.R | \
        sed "s/SELECT_NUMBERS/$SELECT_NUMBERS/g" | \
        sed "s/RATIO/$RATIO/g" | > $work_dir/WGCNA_$RATIO.R

    nohup Rscript $work_dir/WGCNA_$RATIO.R > $work_dir/WGCNA_$RATIO.log &

done
```
**We select the vsd_100 as final network consruction** as shown under the logfile 
```
Allowing parallel execution with up to 30 working processes.
pickSoftThreshold: will use block size 2005.
 pickSoftThreshold: calculating connectivity for given powers...
   ..working on genes 1 through 2005 of 22310
   ..working on genes 2006 through 4010 of 22310
   ..working on genes 4011 through 6015 of 22310
   ..working on genes 6016 through 8020 of 22310
   ..working on genes 8021 through 10025 of 22310
   ..working on genes 10026 through 12030 of 22310
   ..working on genes 12031 through 14035 of 22310
   ..working on genes 14036 through 16040 of 22310
   ..working on genes 16041 through 18045 of 22310
   ..working on genes 18046 through 20050 of 22310
   ..working on genes 20051 through 22055 of 22310
   ..working on genes 22056 through 22310 of 22310
   Power SFT.R.sq  slope truncated.R.sq  mean.k. median.k. max.k.
1      1   0.1350  1.630          0.927 4820.000  4.82e+03 7280.0
2      2   0.0652 -0.688          0.924 1590.000  1.56e+03 3420.0
3      3   0.3490 -1.510          0.958  654.000  6.24e+02 1890.0
4      4   0.5390 -2.070          0.967  308.000  2.83e+02 1150.0
5      5   0.6490 -2.360          0.982  161.000  1.41e+02  746.0
6      6   0.7210 -2.460          0.991   90.800  7.45e+01  508.0
7      7   0.7600 -2.510          0.995   54.500  4.15e+01  359.0
8      8   0.7780 -2.510          0.992   34.400  2.41e+01  261.0
9      9   0.7860 -2.450          0.985   22.800  1.45e+01  194.0
10    10   0.8020 -2.280          0.979   15.700  8.90e+00  148.0
11    12   0.9170 -1.930          0.993    8.170  3.64e+00   99.5
12    14   0.9500 -1.850          0.983    4.760  1.59e+00   82.4
13    16   0.9710 -1.740          0.980    3.030  7.42e-01   70.6
14    18   0.9770 -1.630          0.976    2.080  3.61e-01   61.9
15    20   0.9840 -1.550          0.981    1.510  1.85e-01   56.3
16    22   0.9830 -1.490          0.979    1.150  9.75e-02   52.4
17    24   0.9830 -1.450          0.981    0.911  5.27e-02   49.2
18    26   0.9850 -1.400          0.983    0.742  2.93e-02   46.4
19    28   0.9850 -1.370          0.984    0.620  1.65e-02   43.9
20    30   0.9870 -1.340          0.986    0.528  9.52e-03   41.7
null device 
          1 
 softConnectivity: FYI: connecitivty of genes with less than 16 valid samples will be returned as NA.
 ..calculating connectivities.. 
  scaleFreeRsquared slope
1              0.91 -1.93
null device 
          1 
[1] 12

Allowing parallel execution with up to 30 working processes.
..connectivity..
..matrix multiplication (system BLAS)..
..normalization..
..done.
 ..cutHeight not given, setting it to 0.998  ===>  99% of the (truncated) height range in dendro.
 ..done.
dynamicMods
   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
2646 1156 1156  913  874  861  838  823  805  797  766  612  548  529  515  456 
  16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
 429  425  422  397  374  372  363  356  354  347  336  319  312  311  280  280 
  32   33   34   35   36   37   38   39   40   41   42   43 
 279  261  243  228  205  188  182  163  158  150  141  140 
```
### **Prepare the Formal network construction**
#### Load requried libraries and data input
```r
#install.packages("WGCNA")
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install(c("GO.db", "preprocessCore", "impute"))
install.packages("Rcpp")

### Load libraries
library(WGCNA)
library(cluster)
library(dplyr)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(viridis)
library(matrixStats)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(cowplot)
### emable teh multiple threads
enableWGCNAThreads(nThreads = 30)
```
#### Load expression data and meta-table
```r
#Load expression
TPM <- read.csv("/home/liangyu/1_Project/8_Sorghum/3_WGCNA/Sorghum_20220629_clean_vsd_update96.csv", header = T)
#rename the column name
rownames(TPM) <- TPM$Name
TPM <- TPM[,-1]

#Reorder samples orders by colnames
TPM <- TPM[ , order(names(TPM))]

###Calculate the average numbers for each two replicates
TPM_clean <- data.frame(matrix(ncol = 48, nrow = nrow(TPM)))
vector <- unique(Trait$ID)
rownames(TPM_clean) <- rownames(TPM)
colnames(TPM_clean) <- vector

###Loop calculation for each accesion average
for (i in 1:length(vector)){
  TPM_clean[,grepl(vector[i],colnames(TPM_clean))] <- rowMeans(TPM[,grepl(vector[i],colnames(TPM))])
}
```

#### Load phenotypic and metablomic data
```r
#load oxdative stress table for average
Trait <- read.csv("/home/liangyu/1_Project/8_Sorghum/3_WGCNA/Phenotype/Sorghum_oxidative_stress-Feb2022.csv", 
                  header = T, row.names = 1)
datTrait <- Trait[, 6:34]

datTrait_ave <- data.frame(matrix(ncol = 29, nrow = 48))
vector <- unique(Trait$ID)
colnames(datTrait_ave) <- colnames(datTrait)
rownames(datTrait_ave) <- vector
### Calculate the average by loop
for (i in 1:length(vector)){
  datTrait_ave[grepl(vector[i],rownames(datTrait_ave)),] <- colMeans(datTrait[grepl(vector[i],rownames(datTrait)),])
}
```
```r
#load LI_6800 table for average (morning data)
LI1 <- read.csv("/home/liangyu/1_Project/8_Sorghum/3_WGCNA/Phenotype/LI6800_clean-20220324-morning.csv", 
                  header = T, row.names = 1)
dat_LI1 <- LI1[, 7:18]

LI1_ave <- data.frame(matrix(ncol = 12, nrow = 48))
vector1 <- unique(LI1$ID)
colnames(LI1_ave) <- colnames(dat_LI1 )
rownames(LI1_ave) <- vector1

### Calculate the average by loop
for (i in 1:length(vector1)){
  LI1_ave[grepl(vector1[i],rownames(LI1_ave)),] <- colMeans(dat_LI1[grepl(vector1[i],rownames(dat_LI1)),])
}
```
```r
#load LI_6800 table for average (afternoon data)
LI2 <- read.csv("/home/liangyu/1_Project/8_Sorghum/3_WGCNA/Phenotype/LI6800_clean-20220324-afternoon.csv", 
                  header = T, row.names = 1)
dat_LI2 <- LI2[, 7:18]

LI2_ave <- data.frame(matrix(ncol = 12, nrow = 48))
vector2 <- unique(LI2$ID)
colnames(LI2_ave) <- colnames(dat_LI2)
rownames(LI2_ave) <- vector2

### Calculate the average by loop
for (i in 1:length(vector2)){
  LI2_ave[grepl(vector2[i],rownames(LI2_ave)),] <- colMeans(dat_LI2[grepl(vector2[i],rownames(dat_LI2)),])
}
```
#### Load genes with high MAD score
```r
### Take the top 85% of genes

datExpr <- as.data.frame(t(TPM_clean[order(apply(TPM_clean,1,mad), decreasing = T)[1:19895],]))
#write.csv(data.frame(transcripts_ID = colnames(datExpr)),"/home/liangyu/1_Project/8_Sorghum/3_WGCNA/MAD85_TranscriptsID.csv", quote=F, row.names = F)

### Transform TPM into Z-score
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

Zscore_TPM <- scale_rows(as.data.frame(t(datExpr)))
write.csv(Zscore_TPM, "/home/liangyu/1_Project/8_Sorghum/3_WGCNA/Zscore_TPM_20220331.csv", quote = F)
```

#### Calculate the soft thresholding power (beta)
The cutoff of correlation is 0.8
```r
enableWGCNAThreads(nThreads = 10)

powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

 # Plot the results:

  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  sft$powerEstimate
  
k <- softConnectivity(datE=datExpr,power=sft$powerEstimate) 
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
```

#### Step-by-step network construction and module detection
```r
#Co-expression similarity and adjacency
softPower = sft$powerEstimate;
adjacency = adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);

#set the minimum module size
minModuleSize = 100;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 4, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
```
#### Merging of modules whose expression profiles are very similar
```r
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

### Chose the cutoff for tree height
### Corrssponded to the genes with 0.8 corelation coefficient will be merged into one module

MEDissThres = 0.2

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

### Plot the final modules been merged 
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

### Export the count data
Module_count <- as.data.frame(table(mergedColors))
write.csv(Module_count, "/mnt/Knives/1_data/Sorghum_WGCNA/Results/MAD85/MODULE.COUNT_MAD85.CSV")
```

#### Compare the epigengene values across samples from each module
```r
#Build the epignegene value matrix for plot
MEs0 <- moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs0$Sample <- rownames(MEs0)
MEs0$Name <- rownames(MEs0)
MEs <- separate(MEs0, Name, c("Accession", "Condition", "Time"))

### Plot Eigengene value for each module in a line graph mode
for (i in 1:nrow(Module_count)){
  plot1 <- ggplot(MEs[MEs$Condition =="WL",],aes(x = Time,y = get(paste0("ME", Module_count[i,1])), color = Accession)) +
  geom_point(size = 1) +
  geom_line(alpha = 0.3,size = 1, aes(group= Accession)) +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90)) + ylab(paste0("WL Eigengene value of ", Module_count[i,1]))

  plot2 <- ggplot(MEs[MEs$Condition =="WW",],aes(x = Time,y = get(paste0("ME", Module_count[i,1])), color = Accession)) +
  geom_point(size = 1) +
  geom_line(alpha = 0.3,size = 1, aes(group= Accession)) +
  theme_bw(base_size = 10) +
  theme(axis.text.x=element_text(angle=90)) + ylab(paste0("WW Eigengene value of ", Module_count[i,1]))
  
  ###  combine plot and export into pdf
  pdf(paste0("/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Eigene-value/eigengene_", Module_count[i,1], ".pdf", sep = ""), width = 8, height = 4)
  #print(plot_grid(plot1, plot2, labels = c('A', 'B'), label_size = 12))
  dev.off()
}
```
#### Module-trait membership analysis
```r
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

moduleColors <- mergedColors
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

corType = "pearson"
robustY = ifelse(corType=="pearson",T,F)
```
#### Oxdative stress data correlation with expression
```r
##Build correlation between oxdative status data and expression module
if (corType=="pearson") {
  modOXCor = cor(MEs, datTrait_ave[,colnames(datTrait_ave)[1:29]], use = "p")
  modOXP = corPvalueStudent(modOXCor, nSamples)
} else {
  modOXCorP = bicorAndPvalue(MEs_col, colnames(datTrait_ave[1:29]) , robustY=robustY)
  modOXCor = modOXCorP$bicor
  modOXP   = modOXCorP$p
}

textMatrix1 = paste(signif(modOXCor, 2), "\n(", signif(modOXP, 1), ")", sep = "")
dim(textMatrix1) = dim(modOXCor)

modOXCor_results <- as.data.frame(modOXCor)
modOXP_results <- as.data.frame(modOXP)

min(modOXCor_results)
max(modOXCor_results)
write.csv(modOXCor_results,
          "/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Phenotype/2_Oxdative_vsd100_Correlation.csv", quote = F)

write.csv(modOXP_results,
          "/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Phenotype/2_Oxdative_vsd100_Pvalue.csv", quote = F)
```
#### Build correlation between LI6800 data and expression module (morning data/afternoon data)
```r
###Build correlation between LI6800 data and expression module (morning)

if (corType=="pearson") {
  modLI1Cor = cor(MEs, LI1_ave[,colnames(LI1_ave)[1:12]], use = "p")
  modLI1P = corPvalueStudent(modLI1Cor, nSamples)
} else {
  modLI1CorP = bicorAndPvalue(MEs_col, colnames(LI1_ave[1:12]) , robustY=robustY)
  modLI1Cor = modLI1CorP$bicor
  modLI1P   = modLI1CorP$p
}

textMatrix2 = paste(signif(modLI1Cor, 2), "\n(", signif(modLI1P, 1), ")", sep = "")
dim(textMatrix2) = dim(modLI1Cor)

modLI1Cor_results <- as.data.frame(modLI1Cor)
modLI1P_results <- as.data.frame(modLI1P)
min(modLI1Cor_results)
max(modLI1Cor_results)
write.csv(modLI1Cor_results, "/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Phenotype/2_LIM_vsd100_correlation.csv", quote = F)
write.csv(modLI1P_results, "/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Phenotype/2_LIM_vsd100_Pvalue.csv", quote = F)

###Build correlation between LI6800 data and expression module (afternoon)
if (corType=="pearson") {
  modLI2Cor = cor(MEs, LI2_ave[,colnames(LI2_ave)[1:12]], use = "p")
  modLI2P = corPvalueStudent(modLI2Cor, nSamples)
} else {
  modLI2CorP = bicorAndPvalue(MEs_col, colnames(LI2_ave[1:12]) , robustY=robustY)
  modLI2Cor = modLI2CorP$bicor
  modLI2P  = modLI2CorP$p
}

textMatrix3 = paste(signif(modLI2Cor, 2), "\n(", signif(modLI2P, 1), ")", sep = "")
dim(textMatrix3) = dim(modLI2Cor)

modLI2Cor_results <- as.data.frame(modLI2Cor)
modLI2P_results <- as.data.frame(modLI2P)
min(modLI2Cor_results)
max(modLI2Cor_results)

write.csv(modLI2Cor_results, "/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Phenotype/2_LIA_vsd100_correlation.csv", quote = F)
write.csv(modLI2P_results, "/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Phenotype/2_LIA_vsd100_Pvalue.csv", quote = F)
```

#### Calculate the gene-module membership (MM)
```r
# names (colors) of the modules
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#Export the gene-module membership (MM table includes all modules)
write.csv(geneModuleMembership, "/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Phenotype/CorrelationGene_module_membership_vsd100.csv")
```

#### GS_MM based on particular trait
```r
### extract the GS of oxidative stress phenotype
Oxidative_list <- colnames(datTrait_ave)
Oxidative_GS <- data.frame(matrix(ncol = length(Oxidative_list) , nrow = ncol(datExpr)))

for (i in 1:length(Oxidative_list)){
  metab <- datTrait_ave[,colnames(datTrait_ave) == Oxidative_list[i]]
  Oxidative_GS[,i] <- as.data.frame(cor(datExpr, metab, use = "p"))
  colnames(Oxidative_GS)[i] <- Oxidative_list[i]
  rownames(Oxidative_GS) <- colnames(datExpr)
}

### export phenotype data
write.csv(Oxidative_GS, "/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Phenotype/Oxidative_GS.csv", quote = F)


### extract the GS of LIC6800 phenotype
LIM_list <- colnames(LI1_ave)
LIM_GS <- data.frame(matrix(ncol = length(LIM_list) , nrow = ncol(datExpr)))

for (i in 1:length(LIM_list)){
  metab <- LI1_ave[,colnames(LI1_ave) == LIM_list[i]]
  LIM_GS[,i] <- as.data.frame(cor(datExpr, metab, use = "p"))
  colnames(LIM_GS)[i] <- LIM_list[i]
  rownames(LIM_GS) <- colnames(datExpr)
}

### export phenotype data
write.csv(LIM_GS, "/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Phenotype/LIM_GS.csv", quote = F)
```

#### Export network for cytoscape
```r
for (i in 1:nrow(Module_count)){
  # Select module
  module <- Module_count[i,1]
  # Select module probes
  probes = colnames(datExpr)
  moduleColors = mergedColors
  inModule = (moduleColors==module);
  modProbes = probes[inModule];
  # Select the corresponding Topological Overlap
  modTOM <- TOM[inModule, inModule]
  dimnames(modTOM) <- list(modProbes, modProbes)

   cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile = 
  paste0("/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Module/edges-vsd100", paste(module), ".csv"),
  nodeFile = 
  paste0("/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Module/nodes-vsd100", paste(module), ".csv"),
  weighted = TRUE,
  threshold = 0,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule])
}

save.image("/mnt/Knives/1_data/Sorghum_WGCNA/Sorghum_TOP85MADscore.RData")
```

### **Slection of trait-associated genes**

```r
library(ggplot2)
library(dplyr)
```
#### Load MM and GS files
```r
LIM <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/GS_MM_selection/LIM_GS.csv", header = T)
colnames(LIM)[1] = "id"
Oxidative <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/GS_MM_selection/Oxidative_GS.csv", header = T)
colnames(Oxidative)[1] = "id"

MM <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/GS_MM_selection/MM_MAD85_clean.csv", 
               header = T, row.names = 1)

### lable the module information under the GS files
LIM_clean <- left_join(MM, LIM, by = "id")
Oxidative_clean <- left_join(MM, Oxidative, by = "id")
```

#### Load gene annotation information 
```r
GeneID <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/MAD85_functional_annocation.csv",header = T)
```
#### Load metatable containing the interested module-trait combinations
```r
Meta1 <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/GS_MM_selection/MM_oxidative_meta.csv", header = T)
Meta2 <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/GS_MM_selection/MM_LIC_meta.csv", header = T)


### Create sub gene list 
for (i in 1:nrow(Meta1)){
  temp_df <- Oxidative_clean[Oxidative_clean$module == Meta1[i,1],c("id","MM_score",Meta1[i,2])]
  temp_df <- temp_df[order(temp_df[,3]),]
  ### Filter by GS and MM
  temp_df$ABS_MM <- abs(temp_df[,2])
  temp_df$ABS_GS <- abs(temp_df[,3])
  temp_df$Group <- "FALSE"
  temp_df[temp_df$ABS_MM > 0.7 & temp_df$ABS_GS > 0.5,]$Group = "TRUE"
  temp_clean <- temp_df[temp_df$ABS_MM > 0.7 & temp_df$ABS_GS > 0.5,]
  temp_clean <- left_join(temp_clean,GeneID, by = "id")
  
  write.csv(temp_clean, 
            paste0("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/GS_MM_selection/Results/",Meta1[i,1], "_",Meta1[i,2], ".csv"), quote = F, row.names = F)
  assign(paste0(Meta2[i,1], "_",Meta2[i,2]), temp_df)
}


for (i in 1:nrow(Meta2)){
  temp_df <- LIM_clean[LIM_clean$module == Meta2[i,1],c("id","MM_score",Meta2[i,2])]
  temp_df <- temp_df[order(temp_df[,3]),]
  ### Filter by GS and MM
  temp_df$ABS_MM <- abs(temp_df[,2])
  temp_df$ABS_GS <- abs(temp_df[,3])
  temp_df$Group <- "FALSE"
  temp_df[temp_df$ABS_MM > 0.7 & temp_df$ABS_GS > 0.3,]$Group = "TRUE"
  temp_clean <- temp_df[temp_df$ABS_MM > 0.7 & temp_df$ABS_GS > 0.4,]
  temp_clean <- left_join(temp_clean,GeneID, by = "id")
  write.csv(temp_clean, 
            paste0("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/GS_MM_selection/Results/",Meta2[i,1], "_",Meta2[i,2], ".csv"), quote = F, row.names = F)
  
  assign(paste0(Meta2[i,1], "_",Meta2[i,2]), temp_df)
}
```

#### Plot MM-GS scater plot with cutoff marked
```r
### Manual plot samples
ggplot(data=darkmagenta_H2O_r, aes(x=ABS_MM, y=ABS_GS, color=Group)) +
  geom_point(size=1.5)+
  scale_colour_manual(values=c("grey60", "#DE6757")) + 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  
  labs(x="Module Membership in darkmagenta module", y="Gene significance for H2O-r")+
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +
  theme(legend.position = 'none')+
  geom_hline(aes(yintercept=0.5),colour="#5B9BD5",lwd=1,linetype=5)+
  geom_vline(aes(xintercept=0.7),colour="#5B9BD5",lwd=1,linetype=5)
```
### **Filter the weight-score of edge data and extract nodes**
title: "Filter_weight"\
author: "A, Nelson; L, Yu"\
date: "3/29/2022"\
output: html_document

**list all samples**
```bash
ls edges-vsd100* > sample.list
```
example of the list
```
edges-vsd100black.csv
edges-vsd100blue.csv
edges-vsd100brown.csv
edges-vsd100cyan.csv
edges-vsd100darkgreen.csv
edges-vsd100darkgrey.csv
edges-vsd100darkmagenta.csv
edges-vsd100darkolivegreen.csv
edges-vsd100darkorange.csv
edges-vsd100darkred.csv
edges-vsd100darkturquoise.csv
edges-vsd100green.csv
edges-vsd100greenyellow.csv
```
**filter the top10% of the connections**
```r
sample.list <- read.table("/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Module/sample.list", header = F)

for (i in 1:nrow(sample.list)){
  sample.df <- read.csv(paste0("/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Module/",
                               sample.list[i,1]), header = T, sep = "\t")
  ### sort weight data from largest to smallest
  Edge_sort <- sample.df[order(-sample.df$weight),]
  
  ### Take the 10% top edge connections
  Edge_clean <- unique(Edge_sort[1: nrow(Edge_sort)*0.1,])
  write.csv(Edge_clean, paste0("/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Module/Clean-edge/",sample.list[i,1]), row.names = F, quote = F)
  
  Node1 <- Edge_clean[,1, drop = F]
  Node2 <- Edge_clean[,2, drop = F]
  colnames(Node1) = "Transcripts"
  colnames(Node2) = "Transcripts"
  
  ### Build non-duplciated transcripts list 
  temp <- unique(bind_rows(Node1, Node2))
  Edge_count <- data.frame(table(bind_rows(Node1, Node2)$Transcripts))
  
  Node_clean <- data.frame(Node =  temp$Transcripts, Module = sample.list[i,1])
  Node_clean$Module <- str_replace(Node_clean$Module, "edges-MAD85", "") 
  Node_clean$Module <- str_replace(Node_clean$Module, ".csv", "")
  
  write.csv(Node_clean, paste0("/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Module/Clean_node/",sample.list[i,1]), row.names = F, quote = F)
  write.csv(Edge_count, paste0("/mnt/Knives/1_data/Sorghum_WGCNA/Results/vsd100/Module/Edge_count/",sample.list[i,1]), row.names = F, quote = F)
}
```
#### Match module membership (MM) to filtered nodes
```r
MM_score <- read.csv("/home/liangyu/1_Project/8_Sorghum/3_WGCNA/CorrelationGene_module_membership_MAD85.csv", header= T)
colnames(MM_score)[1] <- "id"
colnames(MM_score) <- str_replace(colnames(MM_score), "MM", "")
Module_color <- colnames(MM_score[,-1])

### Load ids which associated with the module color
TranscriptsID <- read.table("/home/liangyu/1_Project/8_Sorghum/3_WGCNA/Transcripts.matchinglist.txt", header = T)

#Extracting certain module
for (i in 1:length(Module_color)){
  ID_temp <- TranscriptsID[TranscriptsID$module == Module_color[i], ]
  MM_temp <- MM_score[,c("id", Module_color[i])]
  MM_summary <- left_join(ID_temp,MM_temp, by = "id")
  colnames(MM_summary)[3] <- "MM_score"
  
  assign(paste0(Module_color[i], "_MM_summary"),MM_summary)
}

MM_list <- lapply(ls(pattern = "_MM_summary"), get)
MM_combine <- bind_rows(MM_list)

write.csv(MM_combine, "/home/liangyu/1_Project/8_Sorghum/3_WGCNA/MM_MAD85_clean.csv")
```

## STEP 5 GO and KEGG enrichment analysis for different sets of genes
title: "3_GO_Enrichment_Sorghum"\
author: "Li'ang Yu Andrew Nelson"\
date: "2021. 10. 24"\
output: html_document

This pipeline is aiming used for batch processing for GO and KEGG enrichement derived from diffrent set of genes

### **Install and load packages**
```r
devtools::install_github('GuangchuangYu/clusterProfiler')
BiocManager::install("AnnotationHub")

library(dplyr)
library(ggplot2)
library(reshape2)
library(topGO)
library(plyr)
library(tidyr)
library(enrichplot)
library(scales)
library(clusterProfiler)
library(GOplot)
library(stringr)
```
### **Load database**
```r
### Load genome-wide GO list
geneID2GO <- readMappings("D:/3_Genome/3_Gene_info/Sb_geneid2go.map")
geneNames <- names(geneID2GO)

### Load databse for KEGG
term2gene <- read.table("D:/3_Genome/3_Gene_info/Sb_gene2ko.tab",header=T,sep="\t")
term2name <- read.csv("D:/3_Genome/3_Gene_info/At_k2ko_description.csv", header = T)
```
### **Load different set of genes**
**Set 1: Enrichement of genes from respective WGCNA modules**\
```r
#Load sample_list for reading respective nodes
sample.list <- read.table("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/MAD85/Node/node.list", sep ="\t", header = F)

for (i in 1:nrow(sample.list)){
  sample.df <- read.csv(paste0("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/MAD85/Node/",
                               sample.list[i,1]), header = T)
  if(nrow(sample.df) != 0){
    assign(sample.list[i,1], sample.df)
  } 
}

#Create Meta phenotype sheet by combining muliple files
merge_list <- lapply(ls(pattern = "node"), get)
MOdule_merge <- bind_rows(merge_list)
MOdule_merge <- MOdule_merge[,c(1,2)]
colnames(MOdule_merge) <- c("id","module")
write.csv(MOdule_merge, "C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/MAD85/Node/Module_merge.csv", quote = F)

### generate the unique list of module information
module_colors <- data.frame(module = unique(MOdule_merge$module))
module_colors <- na.omit(module_colors)
```
**Note: the isoform should be removed before performing the enrichment**
Header of the module merge table as shown 
| id | module |
| ---- | ------ |
| Sobic.009G250100.1 | bisque4 |
| Sobic.007G213800.1 | bisque4 |
| Sobic.006G263600.1 | bisque4 |

**Set 2: Enrichement of genes been specifically modified across time-point, accessions, and treatment**\
**Set 3: Enrichement of genes been constantly modified across time-point, accessions, and treatment**\
**Set 4: Enrichement of AME targeted downstream genes**
```r
#Load sample_list for 
sample.list <- read.table("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/16_AME/GO_list", sep ="\t", header = F)

for (i in 1:nrow(sample.list)){
  sample.df <- read.table(paste0("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/16_AME/AME_peak/",
                               sample.list[i,1]), header = F)
  if(nrow(sample.df) != 0){
    assign(sample.list[i,1], sample.df)
  } 
}

merge_list <- lapply(ls(pattern = "clean_GO"), get)
MOdule_merge <- bind_rows(merge_list)
MOdule_merge <- MOdule_merge[,c(1,2)]
colnames(MOdule_merge) <- c("module","Gene")

MOdule_clean <- MOdule_merge[,c(2,1)]
```
Header of the module clean table as shown 
| Gene | module |
| ---- | ------ |
| Sobic.001G009400 | GT1_col_a_m1 |
| Sobic.009G156600 | GT1_col_a_m1 |
| Sobic.009G15690 | GT1_col_a_m1 |

### **Wrap gene function information to each module**
```r
### Load gene function
Gene_Function <- read.csv("D:/3_Genome/3_Gene_info/Sorghun_annotation_info.clean.csv", header = T)
Gene_sub <- data.frame(matrix(ncol = 6, nrow = nrow(MOdule_merge)))
colnames(Gene_sub) <- c("transcriptName", "Best.hit.arabi.name", "arabi.symbol", "arabi.defline", "KO", "GO")

for (i in 1:nrow(MOdule_merge)){
        temp.df <- Gene_Function[Gene_Function$transcriptName == MOdule_merge[i,1],
                c("transcriptName", "Best.hit.arabi.name", "arabi.symbol", "arabi.defline", "KO", "GO")]
        if(nrow(temp.df) != 0 ){
           Gene_sub[i,] <- temp.df
        } else {
        Gene_sub[i,] = "NA"
        }
}

Gene_info <- cbind(MOdule_merge, Gene_sub)
### write to pdf
write.table(Gene_info, 
          "C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/7_WGCNA/MAD85/Functional_annotation.txt", 
          row.names = F, quote = F, sep = "\t")
```
### Perform functional enrichment
For three sets of data
```r
### for module enrichement
module_colors <- data.frame(color = unique(MOdule_clean$module)) 

### for AME peack analysis
module_colors <- data.frame(color = unique(MOdule_merge$module)) 

### for Modification analyis
MOdule_clean <- read.csv("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/11_MODs_profile/1_MODs_transcripts_GO-enrichemnt_May10-2022.csv", #header = T)
```
Start the enrichment analysis
```r
module_colors <- data.frame(color = unique(MOdule_clean$module)) 

for (i in 1:nrow(module_colors)){
  myInterestingGenes <- MOdule_clean[MOdule_clean$module == module_colors[i,1],1]
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames

  
  #Build the topGOdata matrix for each ontology categories
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=5)

  GOdata_CC <- new("topGOdata", ontology = "CC", allGenes = geneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO,nodeSize=5)

  GOdata_MF <- new("topGOdata", ontology = "MF", allGenes = geneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=5)
  
  #Perform enrichemnt analysis 
  resultTopGO.BP <- runTest(GOdata_BP, algorithm = "classic", statistic = "fisher")
  resultTopGO.CC <- runTest(GOdata_CC, algorithm = "classic", statistic = "fisher" )
  resultTopGO.MF <- runTest(GOdata_MF, algorithm = "classic", statistic = "fisher" )
  
  allBPGO = genesInTerm(GOdata_BP)
  allCCGO = genesInTerm(GOdata_CC)
  allMFGO = genesInTerm(GOdata_MF)
  
  ### Load function for extracting genes associated with a GO ID
  BP_ANOTATION = lapply(allBPGO,function(x) x[x %in% myInterestingGenes])
  CC_ANOTATION = lapply(allCCGO,function(x) x[x %in% myInterestingGenes])
  MF_ANOTATION = lapply(allMFGO,function(x) x[x %in% myInterestingGenes])
  
  
  #Transform data and save into csv
  BP_Res <- GenTable(GOdata_BP, classicFisher = resultTopGO.BP,  topNodes = 20)
  BP_Res$Class = "Biological_Process"
  for (x in 1:nrow(BP_Res)){
    BP_Res[x,8] <- str_c(BP_ANOTATION[[BP_Res[x,1]]],collapse = ",")
  }
  
  CC_Res <- GenTable(GOdata_CC, classicFisher = resultTopGO.CC,  topNodes =  20)
  CC_Res$Class = "Cellular_component"
  for (x in 1:nrow(CC_Res)){
    CC_Res[x,8] <- str_c(CC_ANOTATION[[CC_Res[x,1]]],collapse = ",")
  }

  MF_Res <- GenTable(GOdata_MF, classicFisher = resultTopGO.MF,  topNodes =  20)
  MF_Res$Class = "Molecular_function"
  for (x in 1:nrow(MF_Res)){
    MF_Res[x,8] <- str_c(MF_ANOTATION[[MF_Res[x,1]]],collapse = ",")
  }

  #Combine three types of GO for output
  GO_combine <- rbind(MF_Res,CC_Res,BP_Res)
  GO_combine$classicFisher <- as.numeric(GO_combine$classicFisher)
  
  ### filter enriched term using adjusted P-value
  GO_combine <- GO_combine[GO_combine$classicFisher < 0.05 ,] 
  names(GO_combine)[8] = "GeneList"

  write.csv(GO_combine, 
            paste0("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/8_GO/MAD85/GO_MAD85_",module_colors[i,1],".csv"), sep ="\t",
            quote = F, row.names = F)
  #write.csv(GO_combine, 
            #paste0("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/16_AME/",module_colors[i,1],"peak_GO.csv"), sep ="\t",
            #quote = F, row.names = F)
   
  #write.csv(GO_combine, 
            #paste0("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/11_MODs_profile/MODs_GO/",module_colors[i,1],"specific_GO.csv"), sep ="\t",
            #quote = F, row.names = F)
}
```

### **Load KEGG pathway data**
```r
##Import gene list
gene <- read.csv("D:/3_Genome/3_Gene_info/Gh_module_gene.csv",header = F,sep=",")
gene <- as.factor(gene$V1)

# import KEGG term gene
K2gene <- read.table("D:/3_Genome/3_Gene_info/Sb_geneid2KEGG.map",header=T,sep="\t")
K2Ko <- read.table("D:/3_Genome/3_Gene_info/Gh_k2ko.map",header=T,sep="\t")

gene2ko=merge(K2gene,K2Ko,by="K")
write.table(gene2ko[,c(3,2)],"D:/3_Genome/3_Gene_info/Sb_gene2ko.tab",row.names = F,sep = "\t")
```

### **Perform the KEGG enrichement analysis**
```r
for (i in 1:nrow(module_colors)){
  KEGG_gene <- as.factor(MOdule_clean[MOdule_clean$module == module_colors[i,1],1])
  KEGG <- enricher(KEGG_gene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, 
               pAdjustMethod = "BH",qvalueCutoff = 0.05)
  write.csv(na.omit(data.frame(KEGG)), 
            paste0("C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/9_KEGG/MAD85/KEGG_MAD85_",module_colors[i,1],".csv"), row.names = F)
}
```

## STEP 6 Enrichment test of binding motifs from promoter regions
The pipeline is used to perform motif enrichment of 1kb promoter sequences based on DAPseq data
[**DAPseq database**](http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php?)

### **Directories information for STEP 6**
```bash
genome_DIR="/home/liangyu/1_Project/8_Sorghum/6_Annotation"
work_dir="/home/liangyu/1_Project/8_Sorghum/7_Motif"
seq_dir="/mnt/Knives/1_data/Sorghum_WGCNA/Results/MAD85/Clean"
```

### **Extract promoter sequences baesd on sorghum annotation**
Here's the head of Gene_model.list: (Only primary trascnripts included)
```
Sobic.007G071500.1
Sobic.008G049200.2
Sobic.001G149500.1
Sobic.003G431700.1
Sobic.004G051100.1
Sobic.007G136900.1
Sobic.002G211000.1
Sobic.002G352100.1
Sobic.010G000500.1
Sobic.002G030200.1
```

```bash
### extract metafeature based on mRNA annotation
grep 'mRNA' $genome_DIR/Sorghum.gff_primary_transcripts.gff3 > $genome_DIR/Sorghum_PT-mRNA.gff3

for ID in $(cat $work_dir/Gene_model.list);
do
	grep "$ID" $genome_DIR/Sorghum_PT-mRNA.gff3 | cut -f 1,4,5,9 >> $work_dir/coordindate.list  

done

cut -d ':' -f1 coordindate.list > coordindate_clean.list 
sed -i 's/ID=//g' coordindate_clean.list 

### Extract the 1kb promoter sequences based on selected transcripts
cat coordindate_clean.list |awk -F ' ' '{print($2-1000)}' > temp
paste coordindate_clean.list temp > coordinate.final
awk '{ print $1 ":" $5 "-" $2}' coordinate.final >  coordinate.txt

### Use samtools to extract sequences
for LINE in $(cat $work_dir/coordinate.txt);
do 
	samtools faidx $genome_DIR/Sbicolor_454.fa $LINE >> Promoter_seq.fasta
done

### Combine the coordinate.txt and transcipts ID into the coordinate list 
paste $work_dir/coordinate.txt $work_dir/Gene_model.list > $work_dir/format.list
```

### **Replace the transcript ID to coordinate list**
Here's the head of format.list
```Chr07:7913921-7914921	Sobic.007G071500.1
Chr08:4847187-4848187	Sobic.008G049200.2
Chr01:12012864-12013864	Sobic.001G149500.1
Chr03:73420119-73421119	Sobic.003G431700.1
Chr04:4107527-4108527	Sobic.004G051100.1
Chr07:56438996-56439996	Sobic.007G136900.1
Chr02:60381768-60382768	Sobic.002G211000.1
Chr02:71555302-71556302	Sobic.002G352100.1
Chr10:30712-31712	Sobic.010G000500.1
```

```bash
IFS=$'\n';

for LINE in $(cat $work_dir/format.list);
do 
	A=$(echo ${LINE} | awk '{ print $1}')
	B=$(echo ${LINE} | awk '{ print $2}')
	sed "s/$A/$B/g" Promoter_seq.fasta > Promoter_seq.clean.fasta
done
```
### **Perform the enrichment analyiss based on AME pipeline**
AME link: https://meme-suite.org/meme/tools
download databased "ArabidopsisDAPv1.meme"

Here's the head of module list
```
Clean-nodeedges-MAD85bisque4.csv
Clean-nodeedges-MAD85black.csv
Clean-nodeedges-MAD85brown.csv
Clean-nodeedges-MAD85brown4.csv
Clean-nodeedges-MAD85cyan.csv
Clean-nodeedges-MAD85darkgreen.csv
Clean-nodeedges-MAD85darkgrey.csv
Clean-nodeedges-MAD85darkmagenta.csv
Clean-nodeedges-MAD85darkolivegreen.csv
Clean-nodeedges-MAD85darkorange.csv
Clean-nodeedges-MAD85darkorange2.csv
Clean-nodeedges-MAD85darkred.csv
Clean-nodeedges-MAD85darkslateblue.csv
Clean-nodeedges-MAD85darkturquoise.csv
Clean-nodeedges-MAD85floralwhite.csv
```
```bash
IFS=$'\n';
for LINE in $(cat $work_dir/Module.list);
do 
  NAME=$(echo $LINE | cut -d "." -f1 | cut -d "-" -f3)
	cat ${seq_dir}/$LINE | cut -d "," -f1 | grep -v 'node' > ${work_dir}/$NAME.list 

  ###extract promoter sequences to differnet moudles
  faSomeRecords ${work_dir}/Promoter_seq.fasta ${work_dir}/$NAME.list ${work_dir}/Sequences/$NAME.promoter.fasta
	
  ###Perform motif enrichment analysis 
	ame --verbose 1 --oc ${work_dir}/Results/$NAME --o ${work_dir}/Results/$NAME \
	--scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 \
	--control --shuffle-- --kmer 2 ${work_dir}/Sequences/$NAME.promoter.fasta $work_dir/ArabidopsisDAPv1.meme

  ###revise name of files
	mv ${work_dir}/Results/$NAME/ame.tsv ${work_dir}/Results/${NAME}/$NAME.tsv
	mv ${work_dir}/Results/$NAME/ame.html ${work_dir}/Results/${NAME}/$NAME.html
	mv ${work_dir}/Results/$NAME/sequences.tsv ${work_dir}/Results/${NAME}/$NAME.sequence

		 
  ###Extract the gene ID with DAPseq signal from database 
	awk -v $NAME '{print $0}' $work_dir/AME-MAD85_TF2DAPseq > ${work_dir}/Results/${NAME}/$NAME.TF2DAPseq

	for i in $(cat ${work_dir}/Results/${NAME}/$NAME.TF2DAPseq | cut -f4);
	do
		grep $i ${work_dir}/Results/${NAME}/$NAME.tsv >> ${work_dir}/Results/${NAME}/$NAME.select
	done

  ###extract target sequences based on enriched motifs and only targets marked with "true positive will be kept"
	for i in $(cat ${work_dir}/Results/${NAME}/$NAME.select | cut -f3);
	do 
		grep $i ${work_dir}/Results/${NAME}/$NAME.sequence | awk '$6 == "tp" {print $0}'  >> ${work_dir}/Results/${NAME}/$NAME.peak
	done

  ###Create the clean file with information for cytoscape downstream analysis
	cut -f2,3 ${work_dir}/Results/${NAME}/$NAME.peak | cut -d "." -f2,3 | awk '$3 = "DNA-binding" {print $0}' > ${work_dir}/Results/${NAME}/$NAME.clean

done
```

# **Section III High-throughput data processing for post-transcriptional modifications (PTMs) using RNA-seq data**
## STEP 1 Identificaiton of PTMs using Tophat2 mapping and HAMR
The github page of [**HAMR**](https://github.com/GregoryLab/HAMR) and The [**Rdata libraries**](http://tesla.pcbi.upenn.edu/hamr/)can be accessed here accordingly.  

### **Container information**
```bash
singularity pull docker://reetututeja/hamr_xi:1.4
singularity pull docker://broadinstitute/gatk3:3.5-0
singularity pull docker://quay.io/biocontainers/samtools:1.10--h2e538c0_3
singularity pull docker://quay.io/biocontainers/picard:2.23.8--0
singularity pull docker://quay.io/biocontainers/tophat:2.1.1--py27_3
singularity pull docker://quay.io/biocontainers/sra-tools:2.10.9--pl526haddd2b5_0
singularity pull docker://quay.io/biocontainers/salmon:1.5.0--h84f40af_0
Singularity pull docker://reetututeja/data_processing_notifications:1.0
```

### **Directories information for STEP 1**
```bash
###For alignment
raw_reads_dir="/mnt/Milly/Work/Sorghum_epitranscriptomics_project/Sorghum_seq_third_replacement_round"
trim_reads_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/0_Trim"
genome_dir="/mnt/Milly/Work/Sorghum_epitranscriptomics_project/phytozome/Sbicolor/v3.1.1"
bam_out="/mnt/Knives/1_data/5_MODs_bam/sorghum/3_Alignment"

###For HAMR
Docker_dir="/home/liangyu/6_Docker"
bam_out="/mnt/Knives/1_data/5_MODs_bam/sorghum/3_Alignment"
hamr_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/4_HAMR"
```

### **Index the reference genome**
**Note**
1. Should perform this step use the __-large-index__ option 
2. Fasta file and the index should be placed under the same folder
3. The name of fasta file should be same as the index prefix
   
```bash
bowtie2-build --large-index $genome_dir/assembly/Sbicolor_454_v3.0.1.fa Sbicolor
```

### **Use trimmonmatic to remove adapters and low-quality reads**
```bash
###Trim the low quality reads (Treat as paired end)
for i in S97 S98 S99 S100; do
 trimmomatic PE -threads 10 \
    $raw_reads_dir/${i}-Nelson98_R1.fastq.gz \
    $raw_reads_dir/${i}-Nelson98_R2.fastq.gz \
    ${trim_reads_dir}/${i}_forward_paired.fq.gz ${trim_reads_dir}/${i}_forward_unpaired.fq.gz \
    ${trim_reads_dir}/${i}_reverse_paired.fq.gz ${trim_reads_dir}/${i}_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/home/liangyu/anaconda3/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:2 LEADING:3 TRAILING:3 MINLEN:36

  ### check the numbers of reads within the unpaired reads file
   gzip -l ${trim_reads_dir}/${i}_forward_unpaired.fq.gz | awk 'NR==2 {print $2}'
   gzip -l ${trim_reads_dir}/${i}_reverse_unpaired.fq.gz | awk 'NR==2 {print $2}'

done
```

### **Check the LIB types before alignment**
Input LIBTYPE A to allow Salmon to identify libtype in alignment-based mode. The output should be a string like ISR which is stranded data and read comes from the reverse strand.

```bash

fastq_dir="/home/liangyu/1_Project/8_Sorghum/1_Salmon/0_data"
Ref_dir="/home/liangyu/1_Project/8_Sorghum/1_Salmon/1_Genome"
sample_dir="/home/liangyu/1_Project/8_Sorghum/1_Salmon/0_data"
output_dir="/home/liangyu/1_Project/8_Sorghum/1_Salmon"


###build index for salmon
salmon index -t $Ref_dir/Sbicolor_454_v3.1.1.cds.fa -i $Ref_dir/Sbicolor_454_v3.1.1.salmon

IFS=$'\n';
for LINE in $(cat ${sample_dir}/sample.list);
do

	fq_file=$(echo ${LINE} | awk '{ print $1}')
	salmon_file=$(echo ${LINE} | awk '{ print $2}')

###perform quantification for each sample
	salmon quant -l A -1 $fastq_dir/${fq_file}_R1.fastq.gz -2 $fastq_dir/${fq_file}_R2.fastq.gz -i $Ref_dir/Sbicolor_454_v3.1.1.salmon -o ${output_dir}/${salmon_file}.salmon.count -p 10
done
```

**Check the libraries type based on salmon output under lib_format_counts.json**
```
{
    "read_files": "( /home/liangyu/1_Project/8_Sorghum/1_Salmon/0_data/12397_11887_135000_HLVC7BGXH_1_ATGATTGA_R1.fastq.gz, /home/liangyu/1_Project/8_Sorghum/1_Salmon/0_data/12397_11887_135000_HLVC7BGXH_1_ATGATTGA_R2.fastq.gz )",
    "expected_format": "ISF",
    "compatible_fragment_ratio": 0.7954681018031248,
    "num_compatible_fragments": 10647506,
    "num_assigned_fragments": 13385208,
    "num_consistent_mappings": 15184916,
    "num_inconsistent_mappings": 3920084,
    "MSF": 0,
    "OSF": 11226,
    "ISF": 15184916,
    "MSR": 0,
    "OSR": 242,
    "ISR": 149621,
    "SF": 3634616,
    "SR": 110349,
    "MU": 0,
    "OU": 0,
    "IU": 0,
    "U": 0
}
```

### **Perform reads alignment by Tophat2**
For paired-end dataset, **Treat as single-end and perform two rounds**\
**Keep the paired and unpaired reads from both forward and reverse strand reads** while processing the mapping process

The library type filled under the argument should match with the salmon tested output, as suggested by [**salmon manual**](https://salmon.readthedocs.io/en/latest/library_type.html)

**See details below**
| TopHat | Paired-end | Single-end |
| ------- | ----------- | ----------- |
| -fr-unstranded |	-l IU |	-l U|
| -fr-firststrand |	-l ISR | -l SR|
| -fr-secondstrand | -l ISF | -l SF|

```bash
raw_reads_dir="/mnt/Milly/Work/Sorghum_epitranscriptomics_project"
trim_reads_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/0_Trim"
genome_dir="/mnt/Milly/Work/Sorghum_epitranscriptomics_project/phytozome/Sbicolor/v3.1.1"
bam_out="/mnt/Knives/1_data/5_MODs_bam/sorghum/3_Alignment"
hamr_dir="/home/liangyu/1_Project/8_Sorghum/4_HAMR"

for i in S97 S98 S99 S100; 
do
    mkdir ${bam_out}/${i}_F
    mkdir ${bam_out}/${i}_R 

    READS_NUMBER1=$(gzip -l ${trim_reads_dir}/${i}_forward_unpaired.fq.gz | awk 'NR==2 {print $2}')
    if (($READS_NUMBER1 != 0)) 
    then  
        ### Process the forward reads
        tophat2 -p 20 --read-mismatches 12 \
          --read-edit-dist 12 --max-multihits 10 \
          --b2-very-sensitive --transcriptome-max-hits 10 \
          --library-type fr-secondstrand \
          --GTF $genome_dir/annotation/Sbicolor_454_v3.1.1.gene_exons.gff3 \
          --no-coverage-search \
          --output-dir ${bam_out}/${i}_F \
          $genome_dir/assembly/bowtie2index/Sbicolor \
          ${trim_reads_dir}/${i}_forward_paired.fq.gz,${trim_reads_dir}/${i}_forward_unpaired.fq.gz
    else
        ### Process the forward reads
        tophat2 -p 20 --read-mismatches 12 \
          --read-edit-dist 12 --max-multihits 10 \
          --b2-very-sensitive --transcriptome-max-hits 10 \
          --library-type fr-secondstrand \
          --GTF $genome_dir/annotation/Sbicolor_454_v3.1.1.gene_exons.gff3 \
          --no-coverage-search \
          --output-dir ${bam_out}/${i}_F \
          $genome_dir/assembly/bowtie2index/Sbicolor \
          ${trim_reads_dir}/${i}_forward_paired.fq.gz
    fi

    READS_NUMBER2=$(gzip -l ${trim_reads_dir}/${i}_reverse_unpaired.fq.gz | awk 'NR==2 {print $2}')
    if (($READS_NUMBER2 != 0)) 
    then
        ### Process the reverse reads
        tophat2 -p 20 --read-mismatches 12 \
          --read-edit-dist 12 --max-multihits 10 \
          --b2-very-sensitive --transcriptome-max-hits 10 \
          --library-type fr-secondstrand \
          --GTF $genome_dir/annotation/Sbicolor_454_v3.1.1.gene_exons.gff3 \
          --no-coverage-search \
          --output-dir ${bam_out}/${i}_R \
          $genome_dir/assembly/bowtie2index/Sbicolor \
          ${trim_reads_dir}/${i}_reverse_paired.fq.gz,${trim_reads_dir}/${i}_reverse_unpaired.fq.gz
    else
        ### Process the reverse reads
        tophat2 -p 20 --read-mismatches 12 \
          --read-edit-dist 12 --max-multihits 10 \
          --b2-very-sensitive --transcriptome-max-hits 10 \
          --library-type fr-secondstrand \
          --GTF $genome_dir/annotation/Sbicolor_454_v3.1.1.gene_exons.gff3 \
          --no-coverage-search \
          --output-dir ${bam_out}/${i}_R \
          $genome_dir/assembly/bowtie2index/Sbicolor \
          ${trim_reads_dir}/${i}_reverse_paired.fq.gz
    fi      

    ### Get unique mapped reads
    #### convert bam to sam
    samtools view -h -o ${bam_out}/${i}_F/${i}_F.sam ${bam_out}/${i}_F/accepted_hits.bam
    samtools view -h -o ${bam_out}/${i}_R/${i}_R.sam ${bam_out}/${i}_R/accepted_hits.bam

    ### get the unique mapping reads
    grep -P '^\@|NH:i:1$' ${bam_out}/${i}_F/${i}_F.sam > ${bam_out}/${i}_F/$i.uniqueF.sam
    grep -P '^\@|NH:i:1$' ${bam_out}/${i}_R/${i}_R.sam > ${bam_out}/${i}_R/$i.uniqueR.sam

    ### convert sam to bam and 
    samtools view -bSh ${bam_out}/${i}_F/$i.uniqueF.sam > ${bam_out}/${i}_F/$i.uniqueF.bam
    samtools view -bSh ${bam_out}/${i}_R/$i.uniqueR.sam > ${bam_out}/${i}_R/$i.uniqueR.bam

    ### sort the bam files
    samtools sort ${bam_out}/${i}_F/$i.uniqueF.bam > ${bam_out}/${i}_F/$i.unique.sortF.bam
    samtools sort ${bam_out}/${i}_R/$i.uniqueR.bam > ${bam_out}/${i}_R/$i.unique.sortR.bam

    ### merge two sorted bam files
    samtools merge ${bam_out}/$i.clean.bam ${bam_out}/${i}_F/$i.unique.sortF.bam ${bam_out}/${i}_R/$i.unique.sortR.bam
done
```

### **Modify the format of bam files**
```bash
for i in S97 S98 S99 S100; 
do
  ##Use picard to add replace read group
  picard \
  AddOrReplaceReadGroups -I ${bam_out}/$i.bam -O ${hamr_dir}/$i.RG.bam -ID $i -LB D4 -PL illumina -PU barcode -SM $i

    ### Reorder bam file
    picard ReorderSam \
       INPUT=$hamr_dir/$i.RG.bam \
       OUTPUT=$hamr_dir/$i.RGO.bam \
       SEQUENCE_DICTIONARY=$genome_dir/assembly/PicardDict/Sbicolor.dict

    ### index bam file
    samtools index ${hamr_dir}/$i.RGO.bam ${hamr_dir}/$i.RGO.bam.bai

    ### Resolve splice alignments using GATK-
    singularity run --cleanenv $Docker_dir/gatk3_3.5-0.sif java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \ 
        -T SplitNCigarReads \
        -R $genome_dir/Sbicolor.fa \
        -I ${hamr_dir}/$i.RGO.bam \
        -o ${hamr_dir}/$i.resolvedalig.bam \
        -U ALLOW_N_CIGAR_READS
done
```

### **Perform PTMs detection by HAMR**
**Note: bam files cannot be accesed under the mnt while using the singularity**\
sample.resolvedalig.bam- input BAM file\
genome.fasta- reference genome FASTA file\
Min_read_qual- 30\
Min_read_cov- 50\
Seq_error_rate- 0.01\
Hypothesis- H4\
Max_p- 1\
Max_fdr- 0.05 \
Refpercent- 0.05\

```bash
### Loop the annotation for 100 samples
for i in $(seq -w 1 100); 
do
    ### Perform the HAMR using container
   singularity run --cleanenv $Docker_dir/hamr_xi_1.4.sif \
       -fe ${hamr_dir}/S$i.resolvedalig.bam \
       $genome_dir/Sbicolor.fa \
       models/euk_trna_mods.Rdata $S{i}_sample $S{i} 30 50 0.01 H4 1 0.05 0.05
done
```

### **Perform PTMs detection by MODtect**
Initiall, the customized bed file of transcripts to be annotated should be generated for [**MODtect**](https://github.com/ktan8/ModTect) \
examples are shown as below:
```
Chr01	7528504	7534112	Sobic.001G098100.1
Chr01	73061860	73062909	Sobic.001G454200.1
Chr01	77084626	77086567	Sobic.001G502000.1
Chr01	75758166	75758839	Sobic.001G487500.1
Chr01	46568162	46572072	Sobic.001G262500.2
Chr01	2685066	2687274	Sobic.001G035800.1
Chr01	64946911	64949156	Sobic.001G359800.1
Chr01	77211221	77213039	Sobic.001G503600.2
Chr01	65552504	65562801	Sobic.001G366800.1
```
Resutls for each samples will be saved with multiple columns saved

```bash
sample_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/2_Tophat_bam"
modtect_dir="/mnt/Milly/Work/kylep/ModTect/ModTect_1_7_5_1.py"
genome_dir="/mnt/Milly/Work/Sorghum_epitranscriptomics_project/phytozome/Sbicolor/v3.1.1/assembly" 
output_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/2_Tophat_bam/modtect"

for i in $(ls $sample_dir/*bam*);
do
python $modtect_dir/ModTect_1_7_5_1.py ${i} \
  $genome_dir/Sbicolor.fa \
  Chr01 1 2 \
  --regionFile $sample_dir/sorghum_genes_for_liang2.bed \
  --readlength 151 \
  --threads 12 \
  --minBaseQual 30 \
  --minDepth 10 \
  --label ${output_dir}/${i}
done
```
**Reformat the output from MOdtect**
```bash
Modect_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/2_Tophat_bam/modtect"
ModComp_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Modtect"

### Create bed files based on annotation information
#ls $Modect_dir/*R.merged.bam.combined.txt > $Modect_dir/Mod_sample.txt

for i in $(cat $Modect_dir/Mod_sample.txt);
do
  cat  ${Modect_dir}/$i | cut -f1,2,3 | awk '{print $1,$2,$3,$3,$1":"$2,"+"}' > ${ModComp_dir}/$i.Modtect.bed
  sed -i "s/ /\t/g" ${ModComp_dir}/$i.Modtect.bed
done

cd ${ModComp_dir}
for i in $(ls *Modtect.bed);
do 
  new_name=$(echo $i | cut -d "_" -f5 | awk '{print "S"$1".Modtect.bed"}')
  mv ${ModComp_dir}/$i ${ModComp_dir}/$new_name
done
```

**Prepare the meta-table of rep and sample ID**
```
S51	S3	A_WL_TP1
S98	S63	A_WL_TP3
S73	S25	A_WL_TP5
S100	S87	A_WL_TP7
S57	S9	A_WW_TP1
S69	S21	A_WW_TP3
S81	S33	A_WW_TP5
S93	S45	A_WW_TP7
S53	S5	B_WL_TP1
S65	S17	B_WL_TP3
S74	S26	B_WL_TP5
S89	S41	B_WL_TP7
S56	S8	B_WW_TP1
S68	S20	B_WW_TP3
S80	S32	B_WW_TP5
S92	S44	B_WW_TP7
S54	S6	C_WL_TP1
S66	S18	C_WL_TP3
S78	S28	C_WL_TP5
```
### **Perform the intersection among two replicates of two methods**
**Results of Modtect**
```bash
work_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Modtect"

### Take the intersection of two reps based on rep-meta table
IFS=$'\n';
for LINE in $(cat $work_dir/Meta_table.txt);
do 
  bed_rep1=$(echo ${LINE} | awk '{ print $1}')
  bed_rep2=$(echo ${LINE} | awk '{ print $2 }')
	sample=$(echo ${LINE} | awk '{ print $3 }')

  ### intersect two reps for each  type of modification
	intersectBed -a ${work_dir}/${bed_rep1}.Modtect.bed \
    -b ${work_dir}/${bed_rep2}.Modtect.bed -wo | cut -f1,2,3,4,5,6 | sort | uniq > ${work_dir}/${sample}.intersect.modtect.clean
done
 ### count numbers of mods
for i in $(ls $work_dir);
do 
  wc -l $i >> count.txt; 
done
```

**Results of HAMR**
```bash
work_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/HAMR"
bed_dir="/home/liangyu/1_Project/3_RNA-mods/0_data/Sorghum"

### Take the intersection of two reps based on rep-meta table
IFS=$'\n';
for LINE in $(cat $work_dir/Meta_table.txt);
do 
  bed_rep1=$(echo ${LINE} | awk '{ print $1}')
  bed_rep2=$(echo ${LINE} | awk '{ print $2 }')
	sample=$(echo ${LINE} | awk '{ print $3 }')

  ### reformat bed file
  cat ${bed_dir}/${bed_rep1}.mods.bed | awk '{print $1,$3,$3+1,$1";"$3+1,$5,$6}' > ${work_dir}/bed/${bed_rep1}.HAMR-format.bed
  cat ${bed_dir}/${bed_rep2}.mods.bed | awk '{print $1,$3,$3+1,$1";"$3+1,$5,$6}' > ${work_dir}/bed/${bed_rep2}.HAMR-format.bed

  sed -i 's/ /\t/g' ${work_dir}/bed/${bed_rep1}.HAMR-format.bed
  sed -i 's/ /\t/g' ${work_dir}/bed/${bed_rep2}.HAMR-format.bed

  ### intersect two reps for each  type of modification
	intersectBed -a ${work_dir}/bed/${bed_rep1}.HAMR-format.bed \
    -b ${work_dir}/bed/${bed_rep2}.HAMR-format.bed -wo > ${work_dir}/${sample}.intersect.HAMR.bed
done
```

### **Annotate genes associated with HAMR modified nucleotides**
```bash
note_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/HAMR"
gffdir="/home/liangyu/1_Project/8_Sorghum/5_MODs"

for sample in $(cat $note_dir/sample.list | cut -d "." -f1);
do
    cat ${note_dir}/$sample.intersect.HAMR.clean.bed | sed 's/|/*/g' | sed 's/;/=/g' | awk '{print $1,$2,$3,$4"="$5,$6="0",$7="+"}' > ${note_dir}/$sample.format.bed
    sed -i 's/ /\t/g' ${note_dir}/$sample.format.bed

    cat ${note_dir}/$sample.format.bed | sort | uniq > ${note_dir}/$sample.nodup.bed

    ###run gffcompare to annotate coordinates of mods
    gffcompare -r ${gffdir}/Sorghum.gff_primary_transcripts.gff3 ${note_dir}/$sample.nodup.bed -p MODs -o ${note_dir}/$sample.anno -T 

    ###clean gffcompare annotation file
    cat ${note_dir}/$sample.anno.tracking | cut -f3,4,5 | sed 's/|/ /g' | sed 's/=/ /g' | awk '{print $2,$3,$4,$5,$6}' > ${note_dir}/$sample.anno.clean
    sed -i 's/q1://g' ${note_dir}/$sample.anno.clean

    ###Keep only the C and X scenarios
    awk '$2=="x" || $2=="c" {print $0}' ${note_dir}/$sample.anno.clean > ${note_dir}/$sample.anno.gene.clean

    rm ${note_dir}/$sample.anno.loci
    rm ${note_dir}/$sample.anno.annotated.gtf
    rm ${note_dir}/$sample.anno
    rm ${note_dir}/$sample.anno.tracking
done

for i in $(ls $note_dir/*gene.clean | cut -d "." -f1);
do
  ### calculate the numbers of transcripts been modified
  cat $i.anno.gene.clean | cut -d " " -f1 | sort | uniq -c | wc -l >> 0_combined_gene.count
  ### summary of modified transcripts
  cat $i.anno.gene.clean | cut -d " " -f1 | sort | uniq -c > $i.HAMR.transcripts.count
done
```

### **Annotate genes associated with Modtect modified nucleotides**
```bash
note_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Modtect"
gffdir="/home/liangyu/1_Project/8_Sorghum/5_MODs"

for sample in $(cat $note_dir/sample.list | cut -d "." -f1);
do
    ### reformat the output from Modtect
    cat ${note_dir}/$sample.intersect.modtect.clean |sed 's/:/=/g' | awk '{print $1,$2,$3,$5"="$4,$5="0",$6="+"}' > ${note_dir}/$sample.format.bed
    sed -i 's/ /\t/g' ${note_dir}/$sample.format.bed

    cat ${note_dir}/$sample.format.bed | sort | uniq > ${note_dir}/$sample.nodup.bed

    ###run gffcompare to annotate coordinates of mods
    gffcompare -r ${gffdir}/Sorghum.gff_primary_transcripts.gff3 ${note_dir}/$sample.nodup.bed -p MODs -o ${note_dir}/$sample.anno -T 

    ###clean gffcompare annotation file
    cat ${note_dir}/$sample.anno.tracking | cut -f3,4,5 | sed 's/|/ /g' | sed 's/=/ /g' | awk '{print $2,$3,$4,$5,$6}' > ${note_dir}/$sample.anno.clean
    sed -i 's/q1://g' ${note_dir}/$sample.anno.clean

    ###Keep only the C and X scenarios
    awk '$2=="x" || $2=="c" {print $0}' ${note_dir}/$sample.anno.clean > ${note_dir}/$sample.anno.gene.clean

    rm ${note_dir}/$sample.anno.loci
    rm ${note_dir}/$sample.anno.annotated.gtf
    rm ${note_dir}/$sample.anno
    rm ${note_dir}/$sample.anno.tracking
done

for i in $(ls $note_dir/*gene.clean | cut -f1);
do
  ### calculate the numbers of transcripts been modified
  cat $i.anno.gene.clean | cut -d " " -f1 | sort | uniq -c | wc -l >> 0_combined_gene.count
  ### summary of modified transcripts
  cat $i.anno.gene.clean | cut -d " " -f1 | sort | uniq -c > $i.Modtect.transcripts.count
done
```
### **Intersetion of resutls of modified sites from HAMR and Modtect**
```bash
Modtect_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Modtect"
HAMR_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/HAMR"
final_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Combined"
note_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Annotation"

for sample in $(cat $final_dir/sample.list);
do
    ###intersect results from two methods
		intersectBed -a ${Modtect_dir}/${sample}.intersect.modtect.clean \
      -b ${HAMR_dir}/${sample}.intersect.HAMR.bed -wo | cut -f1,2,3,4 > ${final_dir}/$sample.merge.bed

    ### extract complete information from HAMR bed file using coordinates from merged bed file
    intersectBed -a ${final_dir}/$sample.merge.bed \
      -b ${HAMR_dir}/${sample}.intersect.HAMR.bed -wo | cut -f1,2,3,4,9,10 > ${final_dir}/$sample.final.bed

    awk '{print $1,$2,$3,$1";"$2,$5}' ${final_dir}/$sample.final.bed > ${note_dir}/$sample.clean.bed
done

for i in $(ls $final_dir/*.bed);do wc -l $i >> count.txt; done
```

### **Annotate genes associated with intersected modified nucleotides**
```bash
note_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Annotation"
data_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Combined"
gffdir="/home/liangyu/1_Project/8_Sorghum/5_MODs"

for sample in $(cat $note_dir/sample.list | cut -d "." -f1);
do
    cat ${data_dir}/$sample.merge.bed | sed 's/|/*/g' | sed 's/;/=/g' | awk '{print $1,$2,$3,$4"="$7"="$8,0,$6}' | sed 's/ /\t/g' > ${note_dir}/$sample.format.bed
    cat ${note_dir}/$sample.format.bed | sort | uniq > ${note_dir}/$sample.nodup.bed


    ###run gffcompare to annotate coordinates of mods
    gffcompare -r ${gffdir}/Sorghum.gff_primary_transcripts.gff3 ${note_dir}/$sample.nodup.bed -p MODs -o ${note_dir}/$sample.anno -T 

    ###clean gffcompare annotation file
    cat ${note_dir}/$sample.anno.tracking | cut -f3,4,5 | sed 's/|/ /g' | sed 's/=/ /g' | awk '{print $2,$3,$4,$5,$6,$7}' > ${note_dir}/$sample.anno.clean
    sed -i 's/q1://g' ${note_dir}/$sample.anno.clean

    ###Keep only the C and X scenarios
    awk '$2=="x" || $2=="c" {print $0}' ${note_dir}/$sample.anno.clean > ${note_dir}/$sample.anno.gene.clean

    rm ${note_dir}/$sample.anno.loci
    rm ${note_dir}/$sample.anno.annotated.gtf
    rm ${note_dir}/$sample.anno
    rm ${note_dir}/$sample.anno.tracking
done

for i in $(ls $note_dir/*gene.clean | cut -d "." -f1);
do
  ### calculate the numbers of transcripts been modified
  cat $i.anno.gene.clean | cut -d " " -f1 | sort | uniq -c | wc -l >> combined_gene.count
  ### summary of modified transcripts
  cat $i.anno.gene.clean | cut -d " " -f1 | sort | uniq -c > $i.intersection.transcripts.count
  ### reformat the summary file
  cut -d " " -f3,2 ${note_dir}/$i.intersection.transcripts.count > $i.intersection.transcripts.clean
done
```

### Subtract the information of mapping depth from HAMR
```bash
HAMR_input="/home/liangyu/1_Project/8_Sorghum/5_MODs/1_reads"
Modtect_input=""
depth_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Depth"

for i in $(cat ${HAMR_input}/sample.list | cut -d "." -f1);
do
  awk '{print $1,$2,$2+1,$1";"$2+1,$4,$9,$10,$9+$10,$9/$9+$10}' ${HAMR_input}/$i.mods.txt | sed '/s/ /\t/g' > ${depth_dir}/HAMR/$i.HAMR.depth
done
```

### Subtract the information of mapping depth 
```bash
HAMR_input="/home/liangyu/1_Project/8_Sorghum/5_MODs/1_reads"
Modect_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/2_Tophat_bam/modtect"
depth_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Depth"
combine_dir="/mnt/Knives/1_data/5_MODs_bam/sorghum/7_Mod_comparison/Combined"

### reformat HAMR non-ref and ref read counts
for i in $(cat ${HAMR_input}/sample.list | cut -d "." -f1);
do
 awk '{print $1,$2,$2+1,$1";"$2+1,$4,$9,$10,$9+$10}' ${HAMR_input}/$i.mods.txt | sed 's/ /\t/g' > ${depth_dir}/HAMR/$i.HAMR.depth
 sed -i '1d' ${depth_dir}/HAMR/$i.HAMR.depth
done

### merge the HAMR mapping ratio
IFS=$'\n';
for LINE in $(cat ${depth_dir}/Meta_table.txt);
do
	bed_rep1=$(echo ${LINE} | awk '{ print $1}')
  bed_rep2=$(echo ${LINE} | awk '{ print $2 }')
	sample=$(echo ${LINE} | awk '{ print $3 }')

	intersectBed -a ${depth_dir}/HAMR/$bed_rep1.HAMR.depth \
    		-b ${depth_dir}/HAMR/$bed_rep2.HAMR.depth -wo | cut -f1,2,3,4,5,6,7,8,14,15,16 > ${depth_dir}/HAMR/$sample.HAMR.depth
done

### reformat Modtect ratio file
for i in $(cat $Modect_dir/Mod_sample.txt);
do
  awk '{print $1,$2,$2+1,$1";"$2+1,$3,$8}' ${Modect_dir}/$i > ${depth_dir}/Modtect/$i
  sed -i "s/ /\t/g" ${depth_dir}/Modtect/$i
done
cd ${depth_dir}/Modtect
for i in $(ls *.txt);
do 
 new_name=$(echo $i | cut -d "_" -f5 | awk '{print "S"$1".Modtect.depth"}')
 mv ${depth_dir}/Modtect/$i ${depth_dir}/Modtect/$new_name
done

### merge the modtect mapping ratio
IFS=$'\n';
for LINE in $(cat ${depth_dir}/Meta_table.txt);
do
	bed_rep1=$(echo ${LINE} | awk '{ print $1}')
	bed_rep2=$(echo ${LINE} | awk '{ print $2 }')
	sample=$(echo ${LINE} | awk '{ print $3 }')

	intersectBed -a ${depth_dir}/Modtect/$bed_rep1.Modtect.bed \
    		-b ${depth_dir}/Modtect/$bed_rep2.Modtect.bed -wo | cut -f1,2,3,4,5,6,7,14 > ${depth_dir}/HAMR/$sample.HAMR.depth
done


IFS=$'\n';
for LINE in $(cat ${depth_dir}/Meta_table.txt);
do
	sample=$(echo ${LINE} | awk '{ print $3 }')

  ### merge the stats data derived from two methods
	intersectBed -a ${depth_dir}/Modtect/$sample.Modtect.depth \
    		-b ${depth_dir}/HAMR/$sample.HAMR.depth -wo | cut -f1,2,3,4,5,6,7,13,14,15,16,17,18 > ${depth_dir}/merge/$sample.merge.depth

   ### merge the stats with type of modifications derived from HAMR  
  intersectBed -a ${depth_dir}/merge/$sample.merge.depth \
    		-b ${combine_dir}/$sample.merge.bed -wo | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,20,21 > ${depth_dir}/merge/$sample.merge.depth.mods
done


```

**Metatable format example**
```
bed	depth
S97.mods.bed	S97.mods.txt	
S98.mods.bed	S98.mods.txt	
S99.mods.bed	S99.mods.txt	
S100.mods.bed	S100.mods.txt	
```
The intersection of bed files and gff files was perfomred using [**gffcompare**](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml), only the **C:contained in exon regions** and **X:contained in opposite strand of exons** scenarios were selected for downstream analysis. 

### **Directories information for STEP 2**
```bash
###Scripts for passing parameters for MODs+pipeline.sh
cmddir="/home/liangyu/1_Project/3_RNA-mods/2_scripts"
BEDdir="/home/liangyu/1_Project/8_Sorghum/5_MODs/0_bed"
READdir="/home/liangyu/1_Project/8_Sorghum/5_MODs/1_reads"

gffdir="/home/liangyu/1_Project/8_Sorghum/5_MODs"
outdir="/home/liangyu/1_Project/8_Sorghum/5_MODs"

###meta_table for parameter settins
meta_table_dir="/home/liangyu/1_Project/8_Sorghum/5_MODs"

###direcoty for assiging and cleaning results
cleandir="/home/liangyu/1_Project/8_Sorghum/5_MODs/0_bed_clean"
mods_reads_stat="/home/liangyu/1_Project/8_Sorghum/5_MODs/1_MODs_ratio"
mods_count="/home/liangyu/1_Project/8_Sorghum/5_MODs/2_MODs_count"
```

### **Perform assignement of modified sites (PTMs) to respective adjacent primary transcripts**
**Two R scripts were required to finish the proces**\

**1_MODs_read.R**

```r
#!/usr/bin/Rscript

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(getopt)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
command = matrix(c(
    'reads', 'r', 1, "character",
    'mods', 'm', 1, "character",
    'output_table', 't', 1, "character",
    'output_meta', 'o', 1, "character"), byrow=TRUE, ncol=4)

args = getopt(command)


##help info-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if (!is.null(args$help) || is.null(args$reads) || is.null(args$mods) || is.null(args$output_table) || is.null(args$output_meta)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Load mapping reads file
MODs_mapping <- read.table(args$reads, header = T)
MODs_mapping$bp <- as.integer(MODs_mapping$bp)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Load mods annotation file
MODs_anno <- read.table(args$mods,header = FALSE)
colnames(MODs_anno) <- c("Gene.ID","Class","chr","bp","MODs_type")


## ----------------------------------------------------------merge ratio and loci information to certain gene annotations-----------------------------------------------------------------
#Join dataset (Annotation and SNP table)
##non-genic part will be removed
Read_meta <- left_join(MODs_anno ,MODs_mapping ,by = c("bp","chr"))


#Calculate ratio for each site
Read_meta$ratio <- Read_meta$nonref/Read_meta$total_reads

Read_meta <- Read_meta[Read_meta$Class != "p",]
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(Read_meta, args$output_table, row.names=FALSE, quote = F) 

## -------------------------------------------------------------Calculate modififed ratio for each gene----------------------------------------------------------------------------------
Read_meta_gene_ratio <- Read_meta[,c("Gene.ID","MODs_type","nonref","ref")]%>%
group_by(Gene.ID,MODs_type) %>%
dplyr::mutate(nonref_sum =sum(nonref))%>%
dplyr::mutate(ref_sum =sum(ref))
    

#calculate the ratio of reads under certain genes
Read_meta_gene_ratio$Gene_ratio <- Read_meta_gene_ratio$nonref_sum/(Read_meta_gene_ratio$nonref_sum + Read_meta_gene_ratio$ref_sum)


#Re-structure mods using Dcast function
Read_meta_gene_ratio_dcast <- dcast(unique(Read_meta_gene_ratio[,c("Gene.ID","MODs_type","Gene_ratio")]), Gene.ID~MODs_type)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(Read_meta_gene_ratio_dcast, args$output_meta, row.names=FALSE, quote = F) 
```

**2_MODs_counts.R**
```r
#!/usr/bin/Rscript
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#load basic packages
library(getopt)
library(tidyr)
library(reshape2)
library(dplyr)
## -----------------------------------------------------------------------------------------------------------------------------------
command = matrix(c(
    'ratio', 'r', 1, "character",
    'output1', 'o1', 1, "character",
    'output2', 'o2', 1, "character"), byrow=TRUE, ncol=4)

args = getopt(command)

##help info---------------------------------------------------------------------------------------------------------------------------
if (!is.null(args$help) || is.null(args$ratio) || is.null(args$output1) || is.null(args$output2)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}
## -----------------------------------------------------------------------------------------------------------------------------------
#Load files
Read_Meta <- read.csv(args$ratio, header = T)
## ---------------------------------------calculate mods numbers (all and sep) under same gene-----------------------------------------

MODs_count_sep <- Read_Meta[,c("Gene.ID","MODs_type")]%>%
  group_by(Gene.ID,MODs_type) %>%
  dplyr::mutate(Counts =  n())

#Re-structure mods
MODs_count_sep_dcast <- dcast(MODs_count_sep,  Gene.ID ~ MODs_type)
#calculate all mods numbers under same gene
MODs_count_all <- Read_Meta[,c("Gene.ID","MODs_type")]%>%
  group_by(Gene.ID) %>%
  dplyr::mutate(Counts =  n())

MODs_count_all <- unique(MODs_count_all[,c(1,3)])

## -----------------------------------------------------------------------------------------------------------------------------------
#save results
write.csv(MODs_count_sep_dcast, args$output1, row.names = F, quote = F) 
write.csv(MODs_count_all, args$output2, row.names = F, quote = F) 
```

```bash
# meta_table for parameter settins
meta_table_dir="/home/liangyu/1_Project/8_Sorghum/5_MODs"

#direcoty for assiging and cleaning results
cleandir="/home/liangyu/1_Project/8_Sorghum/5_MODs/0_bed_clean"
mods_reads_stat="/home/liangyu/1_Project/8_Sorghum/5_MODs/1_MODs_ratio"
mods_count="/home/liangyu/1_Project/8_Sorghum/5_MODs/2_MODs_count"
-p MO
IFS=$'\n';
for LINE in $(cat $meta_table_dir/meta_table.txt | grep -v 'depth');do
    ### define variables
    bed_file=$(echo ${LINE} | awk '{ print $1}')
    reads_depth_file=$(echo ${LINE} | awk '{ print $2 }')

    echo $bed_file
    echo $reads_depth_file

    ###perform analysis for each sample
    Rscript $cmddir/0_MODs_Filter.R --reads ${READdir}/$reads_depth_file --mods ${BEDdir}/$bed_file

    ###subtract the SRR ID
    SRR=$(echo ${bed_file} | cut -d '.' -f 1)

    ###Transform the bed file
    cat ${BEDdir}/$bed_file.clean | sed 's/|/*/g' | sed 's/;/=/g' | awk '{print $1,$2,$3,$4"="$5,$6,$7="0"}'| awk '{print $1,$2,$3,$4,$6,$5}' > ${BEDdir}/$SRR.format.bed
    sed -i 's/ /\t/g' ${BEDdir}/$SRR.format.bed

    ###run gffcompare to annotate coordinates of mods
    gffcompare -r ${gffdir}/Sorghum.gff_primary_transcripts.gff3 ${BEDdir}/$SRR.format.bed -p MODs -o ${cleandir}/$SRR.anno -T 

    ###clean gffcompare annotation file
    cat ${cleandir}/$SRR.anno.tracking | cut -f3,4,5 | sed 's/|/ /g' | sed 's/=/ /g' | awk '{print $2,$3,$4,$5,$6}' > ${cleandir}/$SRR.anno.clean
    sed -i 's/q1://g' ${cleandir}/$SRR.anno.clean

    ###Keep only the C and X scenarios
    awk '$2=="x" || $2=="c" {print $0}' ${cleandir}/$SRR.anno.clean > ${cleandir}/$SRR.anno.gene.clean

    ###Remove inter-mediate files
    rm ${cleandir}/$SRR.anno.loci
    rm ${cleandir}/$SRR.anno.annotated.gtf
    rm ${cleandir}/$SRR.anno
    rm ${cleandir}/$SRR.anno.tracking

    ###Merging mapped reads and ratio of modified reads 
    Rscript $cmddir/1_MODs_read.R \
        --reads ${READdir}/$reads_depth_file.clean \
        --mods ${cleandir}/$SRR.anno.gene.clean \
        --output_table ${mods_reads_stat}/$SRR.reads.ratio.csv \
        --output_meta ${mods_reads_stat}/$SRR.reads.stats.csv
    
    ###Generate counts of mods per transcripts
    Rscript $cmddir/2_MODs_counts.R \
        --ratio ${mods_reads_stat}/$SRR.reads.ratio.csv \
        --output1 ${mods_count}/$SRR.counts.sep.csv \
        --output2 ${mods_count}/$SRR.counts.all.csv
done
```

## STEP 3 Characterization of distribution of PTMs across diffrent genomic context and plot the meta-plot
### **Directories information for STEP 3**
```bash
work_dir="/home/liangyu/1_Project/8_Sorghum/8_Metaplot"
bed_dir="/home/liangyu/1_Project/8_Sorghum/8_Metaplot/PTMs"
info_dir="/home/liangyu/1_Project/8_Sorghum/8_Metaplot/Sorghum-info"
```
### **Generate a bed file used for parsing genomic-element feature (gff2bed.R)**
```r
---
title: "gff2bed"
author: "A, Nelson; L, Yu"
date: "4/6/2022"
output: html_document
---
###Load clean file
mRNA <- read.table("/home/liangyu/1_Project/8_Sorghum/8_Metaplot/Sorghum_bed/Sbicolor_mRNA.table")
exon <- read.table("/home/liangyu/1_Project/8_Sorghum/8_Metaplot/Sorghum_bed/Sbicolor_exon.table")
cds <- read.table("/home/liangyu/1_Project/8_Sorghum/8_Metaplot/Sorghum_bed/Sbicolor_cds.table")


###Load gene annoation file
for (i in 1:nrow(mRNA)){
  temp_cds <- cds[cds$V9 == mRNA[i,1],]
  temp_exon <- exon[exon$V9 == mRNA[i,1],]
  mRNA[i,6] <- min(temp_cds$V4)
  mRNA[i,7] <- max(temp_cds$V5)
  mRNA[i,8] <- nrow(temp_exon)
  mRNA[i,9] <- str_c(sort(temp_exon$V4),collapse = ",")
  mRNA[i,10] <- str_c(sort(temp_exon$V5),collapse = ",")
}

mRNA$V11 = paste0(mRNA$V9, ",")
mRNA$V12 = paste0(mRNA$V10, ",")
mRNA$V13 = paste0(mRNA$V1, ".p")

mRNA_clean <- mRNA[,c(1,2,3,4,5,6,7,8,11,12,13)]
colnames(mRNA_clean) <- c("#name", "chrom",	"strand", "txStart", "txEnd", 
                          "cdsStart",	"cdsEnd",	"exonCount","exonStarts","exonEnds","proteinID")

mRNA_clean <- mRNA_clean[order(mRNA_clean$chrom, mRNA_clean$txStart),]

###Export clean files
write.table(mRNA_clean, "/home/liangyu/1_Project/8_Sorghum/8_Metaplot/Sorghum_bed/Sorghum_pred", quote = F, sep = "\t", row.names = F)
write.table(mRNA_clean[mRNA_clean$chrom == "Chr01",], "/home/liangyu/1_Project/8_Sorghum/8_Metaplot/Sorghum_test/Sorghum_Chr01", quote = F, sep = "\t", row.names = F)
```

### **Genrate metafeature count for each sample**
**The head of Metatable for parsing sample information**\
column1: rep1\
column2: rep2\
column3: unified sample name

```bash
S51	S3	A_WL_TP1
S98	S63	A_WL_TP3
S73	S25	A_WL_TP5
S100	S87	A_WL_TP7
S57	S9	A_WW_TP1
S69	S21	A_WW_TP3
S81	S33	A_WW_TP5
S93	S45	A_WW_TP7
S53	S5	B_WL_TP1
S65	S17	B_WL_TP3
S74	S26	B_WL_TP5
S89	S41	B_WL_TP7
S56	S8	B_WW_TP1
S68	S20	B_WW_TP3
S80	S32	B_WW_TP5
S92	S44	B_WW_TP7
S54	S6	C_WL_TP1
S66	S18	C_WL_TP3
S78	S28	C_WL_TP5
```

**Process 96 sorghum samples**
```bash
### Create a clean directory to save merge results from two samples
if [ -d "${bed_dir}/clean" ]
    then
    echo ""
    echo "Directory ${bed_dir}/clean exists. "
else
    mkdir ${bed_dir}/clean exists
fi

### Parse two reps of each samples to for annotation process
IFS=$'\n';
for LINE in $(cat $work_dir/Meta_table.txt);
do 
    bed_rep1=$(echo ${LINE} | awk '{ print $1}')
    bed_rep2=$(echo ${LINE} | awk '{ print $2 }')
	  sample=$(echo ${LINE} | awk '{ print $3 }')

    ### Generate files for each type of PTMs
	for i in m3C t6A Y D m1A m1G m2G;
	do
		grep $i ${bed_dir}/$bed_rep1.mods.bed.clean > ${bed_dir}/${bed_rep1}.${i}.bed
		grep $i ${bed_dir}/$bed_rep2.mods.bed.clean > ${bed_dir}/${bed_rep2}.${i}.bed

        ###intersect two reps for each  type of modification
		intersectBed -a ${bed_dir}/${bed_rep1}.${i}.bed -b ${bed_dir}/${bed_rep2}.${i}.bed -sorted -wo > ${bed_dir}/clean/$sample.$i.merge

        ###Subtract one sets of PTMs from based on merged bed file
		cut -f1,2,3,4,5,6 ${bed_dir}/clean/$sample.$i.merge > ${bed_dir}/clean/$sample.$i.bed

        ###annotates the user supplied BED file containing single nucleotide genomic coordinates of sites of interest.
		intersectBed -a ${bed_dir}/clean/$sample.$i.bed -b $info_dir/Sorghum_sort_pred_annotation -sorted -wo > ${work_dir}/Results/$sample.$i.anno

        ##identifies the region of the transcript in which the user supplied sites fall and converts the transcriptomic coordinates to metagene coordinates
		perl /home/liangyu/bin/rel_and_abs_dist_calc.pl --bed ${work_dir}/Results/$sample.$i.anno --regions $info_dir/region_sizes.txt > ${work_dir}/Results/$sample.$i.dist.measures.txt
	done
done
```

### **Plot metaPlot of each sample or merged samples**
```r
---
title: "Metaplot"
output: html_notebook
---

library("scales")
library("ggplot2")
library("tidyr")

### Diplay color panel
display.brewer.pal(n = 8, name = 'PuOr')
brewer.pal(n = 8, name = "PuOr")

display.brewer.pal(n = 7, name = 'Set2')
brewer.pal(n = 7,name = "Set2")

### Load sample list
sample.list <- read.table("/home/liangyu/1_Project/8_Sorghum/8_Metaplot/Results/measurement/sample.list")

for (i in 1:nrow(sample.list)){
  dist <- read.table(paste0("/home/liangyu/1_Project/8_Sorghum/8_Metaplot/Results/measurement/",
                                   sample.list[i,1],".dist.measures.txt"), header =  T)
  dist$temp <- sample.list[i,1]
  dist <- dist %>% separate(temp, c("Accession","Condition","Time", "PTMs"))
  assign(paste0(sample.list[i,1],"_PTMs"),dist)
}

merge_list <- lapply(ls(pattern = "_PTMs"), get)
PTMs_combine <- bind_rows(merge_list)
PTMs_combine$Sample <- with(PTMs_combine, paste0(Accession, "_",Condition, "_",Time))

### Define class of TPMs to three regions based on rel_location
PTMs_combine$Class <- "fill"
PTMs_combine[PTMs_combine$rel_location < 1,c("Class")] = "5'UTR"
PTMs_combine[PTMs_combine$rel_location >= 1 & PTMs_combine$rel_location < 2,c("Class")] = "CDS"
PTMs_combine[PTMs_combine$rel_location >= 2 & PTMs_combine$rel_location <= 3,c("Class")] = "3'UTR"


### Plot the global pattern of PTMs
ggplot(PTMs_combine[PTMs_combine$PTMs != "m2G",], aes(x=rel_location, color = PTMs)) + 
    labs(x = "Releative position", y = "Density of PTMs") +
    theme_bw() +
    geom_vline(xintercept = 1:2, col = "red",linetype="dotted") +
    geom_freqpoly(binwidth = 0.1) +
    scale_colour_manual(values=brewer.pal(n = 7,name = "Set2")) +
    scale_y_continuous(breaks=seq(0,18000,by=1000))


ggplot(PTMs_combine[PTMs_combine$PTMs != "m2G",], aes(x=utr5_end, color = PTMs)) +  
   labs(x = "position to 5'UTR (bp)", y = "Density of PTMs") +
    xlim(-500,500) + theme_bw() +
    geom_vline(xintercept = 1:2, col = "red",linetype="dotted") +
    geom_freqpoly(binwidth = 30) +
    scale_colour_manual(values=brewer.pal(n = 7,name = "Set2"))

ggplot(PTMs_combine[PTMs_combine$PTMs != "m2G",], aes(x=utr3_st, color = PTMs)) + 
    labs(x = "position to 3'UTR (bp)", y = "Density of PTMs") +
    xlim(-500,500) + theme_bw() +
    geom_vline(xintercept = 1:2, col = "red", linetype="dotted") +
    geom_freqpoly(binwidth = 30) +
    scale_colour_manual(values=brewer.pal(n = 7,name = "Set2"))

### Plot distribution per sample

SampleID <- unique(PTMs_combine[,c("Accession","Condition","Time","Sample")])

for (i in 1:nrow(SampleID)){
 P3 <- ggplot(PTMs_combine[PTMs_combine$Sample == SampleID[i,4],], aes(x=utr3_st, color = PTMs)) + 
    labs(x = "position to 3'UTR (bp)", y = "Density of PTMs") +
    xlim(-500,500) + theme_bw() +
    geom_vline(xintercept = 1:2, col = "red", linetype="dotted") +
    geom_freqpoly(binwidth = 50) +
    ggtitle(SampleID[i,4]) +
    scale_colour_manual(values=brewer.pal(n = 7,name = "Set2"))


  P2 <- ggplot(PTMs_combine[PTMs_combine$Sample == SampleID[i,4],], aes(x=utr5_end, color = PTMs)) +  
   labs(x = "position to 5'UTR (bp)", y = "Density of PTMs") +
    xlim(-500,500) + theme_bw() +
    geom_vline(xintercept = 1:2, col = "red",linetype="dotted") +
    geom_freqpoly(binwidth = 50) +
    ggtitle(SampleID[i,4]) +
    scale_colour_manual(values=brewer.pal(n = 7,name = "Set2"))

  P1 <- ggplot(PTMs_combine[PTMs_combine$Sample == SampleID[i,4],], aes(x=rel_location, color = PTMs)) + 
    labs(x = "Relative position", y = "Density of PTMs") +
    theme_bw() +
    geom_vline(xintercept = 1:2, col = "red",linetype="dotted") +
    geom_freqpoly(binwidth = 0.1) +
    ggtitle(SampleID[i,4]) +
    scale_colour_manual(values=brewer.pal(n = 7,name = "Set2"))
  
  temp_plot <- PTMs_combine[PTMs_combine$Sample == SampleID[i,4],]
  temp_plot$Class <- factor(temp_plot$Class,levels = c("5'UTR", "CDS", "3'UTR"))
  P4 <- ggplot(temp_plot, aes(x=Class, fill = PTMs)) + 
    geom_bar(position = "fill", color = "black") +
    scale_fill_manual(values=brewer.pal(n = 7,name = "Set2")) +
    theme_bw() + 
    ylab("Portions of PTMs")
  
  print(plot_grid(P1, P2, P3, P4, labels = c('A', 'B', 'C', 'D'), label_size = 12))
}

### Compare numbers of PTMs across regions based on condition
Accession_list <- c("A","B","C","D","E","F")

ggplot(PTMs_combine, aes(x=rel_location, color = PTMs)) + 
    labs(x = "Releative position", y = "Density of PTMs") +
    theme_bw() +
    geom_vline(xintercept = 1:2, col = "red",linetype="dotted") +
    geom_freqpoly(binwidth = 0.1) +
    scale_colour_manual(values=brewer.pal(n = 7,name = "Set2")) +
  facet_wrap(~Time)

ggplot(PTMs_combine, aes(x=rel_location, color = PTMs)) + 
    labs(x = "Releative position", y = "Density of PTMs") +
    theme_bw() +
    geom_vline(xintercept = 1:2, col = "red",linetype="dotted") +
    geom_freqpoly(binwidth = 0.1) +
    scale_colour_manual(values=brewer.pal(n = 7,name = "Set2")) +
  facet_wrap(~Condition)

ggplot(PTMs_combine, aes(x=rel_location, color = PTMs)) + 
    labs(x = "Releative position", y = "Density of PTMs") +
    theme_bw() +
    geom_vline(xintercept = 1:2, col = "red",linetype="dotted") +
    geom_freqpoly(binwidth = 0.1) +
    scale_colour_manual(values=brewer.pal(n = 7,name = "Set2")) +
  facet_wrap(~Accession)


ggplot(PTMs_combine[PTMs_combine$Accession == "A",], aes(x=Class, fill = PTMs)) + 
    geom_boxplot(position = "fill", color = "black") +
    scale_fill_manual(values=brewer.pal(n = 7,name = "Set2")) +
    theme_bw() + 
    ylab("Portions of PTMs") +
   facet_wrap(Condition~Time)
```

## STEP 4 Qunatificaiton of modified genes and modified sites across transcriptome
This part of the analysis will be performed using house-hold R scripts
```r
---
title: "Modification_Profile"
author: "Nelson A, Yu L"
date: "2022/2/7"
output: html_document
---

### Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggVennDiagram)
library(ggvenn)
library(ggpubr)
library(UpSetR)
```

### **Reformat the structure of PTMs count**

Head of the MODs_list
```
S3	A_WL_TP1_1			
S51	A_WL_TP1_2			
S63	A_WL_TP3_2			
S98	A_WL_TP3_3			
S25	A_WL_TP5_1			
S73	A_WL_TP5_2			
S87	A_WL_TP7_2			
S100	A_WL_TP7_3			
S9	A_WW_TP1_1			
S57	A_WW_TP1_2
```
Head of the Gene_frame
```
Gene.ID
<chr>
Sobic.004G229300.1				
Sobic.002G064700.1				
Sobic.009G250100.1				
Sobic.009G069800.1				
Sobic.001G259700.1				
Sobic.007G213800.1				
Sobic.001G063600.1				
Sobic.003G050350.1				
Sobic.006G263600.1				
Sobic.001G489600.1
```

```r
### Load modification profiles for comparisons
PATH_MODs ="C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/1_Sample/"
MODs_counts_list = read.table(paste0(PATH_MODs,"MODs_list"), header = F)
Gene_frame <-  read.table(paste0(PATH_MODs,"WGCNA_Transcripts.txt"), header = F)
colnames(Gene_frame) <- "Gene.ID"

for (i in 1:nrow(MODs_counts_list)){
  sample.df <- read.csv(paste0(PATH_MODs,MODs_counts_list[i,1],".counts.all.csv"), header = T)
  MODs.df <- left_join(Gene_frame,sample.df,by = "Gene.ID")
  MODs.df <- as.data.frame(MODs.df[,2])
  colnames(MODs.df)<- MODs_counts_list[i,2]
  
  assign(paste(MODs_counts_list[i,2], "_signature"),MODs.df) 
  assign(paste0(MODs_counts_list[i,1], "_counts"), sample.df)
}

#Combined all DFs to generate meta tables containing levels of modifications
merge_list1 <- lapply(ls(pattern = "signature"), get)
MODs_combine <- bind_cols(merge_list1)
rownames(MODs_combine) <- Gene_frame$Gene.ID 


### subtract meta table of modified genes (genes been modified or not)
#Replace all values not equal with 1
MODs_combine[is.na(MODs_combine)] <- 0


## Use 1 and 0 to reprents the presence and absence of mods
MODs_combine_upset <- MODs_combine
for (i in 1:ncol(MODs_combine_upset)){
MODs_combine_upset[,i] <-  replace(x = MODs_combine_upset[,i], 
                   MODs_combine_upset[,i] != 0.000000	, 
                   values =  1)
}
```
### **Classify total numbers of mods and total genes been modifed**
TM: total modified sits\
TG: total modified genes

```r
### Characterize total numbers of modified sites
MODs_combine[is.na(MODs_combine)] <- 0
TM <- data.frame(TotalMOds = colSums(MODs_combine))
TM$Sample = rownames(TM)
TM$Test = rownames(TM)

TM <- TM %>% separate(Test, c("Accession", "Treatment", "Time", "Rep"))
min(TM$TotalMOds)
max(TM$TotalMOds)

### Characterize total numbers of modified genes
TG<- data.frame(TotalMOds = colSums(MODs_combine_upset))
TG$Sample = rownames(TG)
TG$Test = rownames(TG)

TG <- TG %>% separate(Test, c("Accession", "Treatment", "Time", "Rep"))
TG$Name <- paste0(TG$Accession,"_",TG$Treatment,"_",TG$Time) 

min(TG$TotalMOds)
max(TG$TotalMOds)
```

### **Plot overall patterns of modifed sites across samples**
```r
### Plot TM over time point
my_comparisons <- list(c("TP5","TP1"), c("TP5","TP3"), c("TP5","TP7"))
ggplot(data = TM, aes(x = Time, y = TotalMOds, outlier.size = 0.01, fill = factor(Treatment))) +
  geom_boxplot() + theme_bw() +
  #add jitters 
  geom_jitter(color="red", size=0.4, alpha=0.9) +
  #add custom color 
  #add meam summary information
  stat_summary(fun = mean, color = "darkred", position = position_dodge(0.75),
             geom = "point", shape = 4, size = 3,
             show.legend = FALSE) +
  # add comparisons among three groups
  stat_compare_means(comparisons = my_comparisons,test=t.test, 
                     label.y = c(max(TM$TotalMOds +4000),max(TM$TotalMOds + 1000),max(TM$TotalMOds + 3000))) +
  #add comparisons within each group 
  stat_compare_means(aes(group = Treatment), method = "t.test", label = "p.format") +
  ylab("Total modified sites")


### Plot TM over accessions 
my_comparisons2 <- list(c("A","B"), c("A","D"))

ggplot(data = TM, aes(x = Accession, y = TotalMOds, outlier.size = 0.01, fill = factor(Treatment))) +
  geom_boxplot() + theme_bw() +
  #add jitters 
  geom_jitter(color="red", size=0.4, alpha=0.9) +
  #add custom color 
  #add meam summary information
  stat_summary(fun = mean, color = "darkred", position = position_dodge(0.75),
             geom = "point", shape = 4, size = 3,
             show.legend = FALSE) +
  # add comparisons among three groups
  stat_compare_means(comparisons = my_comparisons2,test=t.test, 
                     label.y = c(max(TM$TotalMOds +1000),max(TM$TotalMOds + 3000))) +
  #add comparisons within each group 
  stat_compare_means(aes(group = Treatment), method = "t.test", label = "p.format") +
  ylab("Total modified sites")
```

### **Characterize the genes been modified genes from two relicas**
0: no modifications\
1: only one rep modified\
2: both reps modified

**Same genes modified from both rep are not restricted to the same sites been modified from both two replicates**

```r
### Calculate the overlapping ratio between replicates among 48 samples 
Sample_name <- unique(TM$Name)
Stats_meta <- data.frame(Name = Sample_name, Rep1 = "fill", 
                         Rep2 = "fill", Overlapped_count = "fill")


### Export the datasheet in terms of presense or absense of mods from two reps
for (i in 1: length(Sample_name)){
  df <- MODs_gene_summary[,grepl(Sample_name[i],colnames(MODs_site_summary))]
  df$Count <- rowSums(df[,c(1,2)])
  colnames(df)[3] <- Sample_name[i]
  
  Stats_meta[i,2] = nrow(df[df[,1]==1,]) 
  Stats_meta[i,3] = nrow(df[df[,2]==1,]) 
  Stats_meta[i,4] = nrow(df[df[,3]==2,]) 
  
  
  ### subtrat the column contanting 0,1,2 note for two reps per sample
  df_rep <- df[,3, drop = F]
  assign(paste0(Sample_name[i], "_SampleRep"), df_rep)
}

### Combine 48 samples with 0,1,2 notation 
Rep_mergeList <- lapply(ls(pattern = "_SampleRep"), get)
Sample_mods_meta <- bind_cols(Rep_mergeList)


write.csv(Sample_mods_meta, "C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/11_MODs_profile/MODs_replica_meta.csv", quote = F)
write.csv(Stats_meta, "C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/11_MODs_profile/MODs_replica.csv", row.names = F, quote = F)

```

### **Export numbers of PTMs  per gene based on averaging two reps**
```r
### Export the datasheet in terms of average mods derived from two reps
for (i in 1: length(Sample_name)){
  df <- MODs_combine[,grepl(Sample_name[i],colnames(MODs_combine))]
  
  ### Conduct if else loop to calculate the mean of mods or assign 
  for (x in 1:nrow(df)){
    if(df[x,1] != 0 & df[x,2] != 0){
      df[x,3] <- rowMeans(df[x,1:2])
      } else {
      df[x,3] <- 0
      }
  } 
  
   colnames(df)[3] <- Sample_name[i] 
  ### subtrat the column contanting 0,1,2 note for two reps per sample
  df_rep <- df[,3, drop = F]
  assign(paste0(Sample_name[i], "_SampleRep_total"), df_rep)
}

Rep_mergeList2 <- lapply(ls(pattern = "_SampleRep_total"), get)
Sample_mods_meta_Total <- bind_cols(Rep_mergeList2)


### Normalize into Z-score based on modified genes
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}


Zscore_mods <- scale_rows(Sample_mods_meta_Total)

write.csv(Sample_mods_meta_Total, "C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/11_MODs_profile/MODs_replica_meta_Total.csv", quote = F)
write.csv(Zscore_mods, "C:/Users/Leon/OneDrive - Cornell University/Projects/9_Sorghum/11_MODs_profile/MODs_replica_meta_Zscore.csv", quote = F)
```
