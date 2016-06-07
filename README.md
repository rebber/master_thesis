# master_thesis
The scripts I used for my master thesis about DNA copy number analysis. 

##File description
This is a brief description of the scripts in this folder, which I wrote and used for my master thesis, and in which order to use them. The "pipeline" presumes that BAM files from targeted sequencing (panel data) are available. If comparison to WGS data is desired, BAM files from this are also necessary. Scripts are named after the main task they are performing, sometimes followed by the data type they are used for and/or some more indication to separate different variants. Many of the scripts are variants of the same process, but adapted for different data types. These have typically been further developed for each new data type I have used them for. Thus, later scripts will be more advanced than earlier versions for the same task. The analysis order of the data types was: Breast v7 (sometimes called only “breast”), MSI-PILOT1 (sometimes MSI1 or in the earliest scripts only MSI), breast v7 PTEN\_del, MSI-PILOT2 (or MSI2), Breast v6, and finally HOPKINS. The scripts will therefore commonly hold more refined and advanced code the later data type they consider. Some scripts originally contained sample IDs (corresponding to patients) inserted into the code. These are confidential due to Swedish law, and therefore in the files here exchanged for toy IDs. Such cases are indicated within the code.

###Run CNVkit
The first thing to do is to run [CNVkit](https://github.com/etal/cnvkit) on the panel data. This is done with the scripts **runCNVkit\*.R**. (Earliest versions are named **CNVkit\*.txt** and contains copies of the code run from the command line.) For this, BAM files from tumor samples are necessary. Some scripts uses normal (or tumor) samples BAM files for reference construction, while other scripts uses already created reference bed files. If a reference is to be constructed a bed file with the panel targets is needed. 

###Compare to WGS
If a comparison with WGS data shall be performed, such data processed with [QDNAseq](http://bioconductor.org/packages/release/bioc/html/QDNAseq.html) is necessary. If QDNAseq files are not available, scripts **runQDNAseq\*.R** are used to create such. QDNAseq data is analyzed in scripts **plotQDNAseq\*.R**.

Comparisons of copy number profiles from panel and WGS data are done with scripts **plotCompWgsPanel\*.R**. 

Gene-wise comparison of CN ratio can also be performed. The first step of this is to calculate the mean segmented log2 CN ratio of each gene from respective data type with scripts **CNVkitSegments2genes\*.R** and **QDNAseqSegments2genes\*.R** (early versions \*.txt). This utilizes the external scripts [qdnaseq2bed.py](https://github.com/dakl/autoseq-scripts/blob/master/qdnaseq2bed.py), [cnvgtf2bed.py](https://github.com/dakl/autoseq-scripts/blob/master/cnvgtf2bed.py) and program [bedtools](https://github.com/arq5x/bedtools2). The comparisons are done through scatter plots created with **scatterPanelVSwgsGene\*.R**.  

###Reference comparison
A set of samples to test different references on is chosen randomly with **saveTestset\*.R**. CNVkit with different references is run on these samples with various **runCNVkit\*.R** scripts, according to naming. 

In order to compare noise levels with MAPD the read count per bin is needed. This is counted with scripts **bedtoolsCountReads\*.R** and then merged with cnr files data through **addReadCounts2Cnr\*.R**. Versions when this is performed within the same script, **bedtoolsCountReadsAdd2Cnr\*.R**, also exist.

CN profile comparisons for different references are done with **plotCompPanelRefs\*.R**. Some of these scripts also save mrpb (median reads per bin) and MAPD data for use in the MAPD comparison. In other cases scripts **saveMAPDmats\*.R** do this. 

####Principal component analysis (PCA)
PCA can be performed to find tumor samples of similar noise structure and correct the noise against these only. In order to do this “smoothed” log2 CN ratio values are needed for each bin. These are calculated with script **addSmooth2cnrRcComb.R**, however this function is in some cases included already in the **bedtoolsCountReadsAdd2Cnr\*.R** scripts. 

The PCA itself and analysis of the results of this is done in scripts **plotPCA\*.R**.  Here clustering of samples and calculation of the median error for each sample as the expected error is performed. Also fitting of a linear function for the expected error based on principal components positions is performed. Subtraction of the expected errors from original log2 CN ratios in the bins is performed to achieve corrected (in the scripts called adjusted) values. **plotPCA\_Bv6AllT\_report.R** saves the PCA figures I used for the report (and also some more figures). Correlations between principal components (and thus also clusters) and parameters from library preparation etc. are analyzed in **analyzeCorClusters.R**.

Script **adjustCnr.R** creates new cnr files with the corrected bin values. **runCNVkitBv6testset\_segmAdjBins.R** uses CNVkit to do segmentation on the corrected cnr files. **plotCompPCAAdjNonadjBv6.R** plots comparisons of copy number profiles from corrected (adjusted) and non-corrected (non-adjusted) log2 CN ratios. 

MAPD and mrpb data for PCA corrected test samples is saved in some cases by **plotCompPCAAdjNonadjBv6.R** and in some cases by **saveMAPDmats\*.R**. 

####Comparison of noise levels using MAPD
The final comparison of remaining noise levels after treatment with the different reference types is done in scripts **plotMAPDcomp\*.R**. Scripts with extension “\_lines” connects the test samples with lines, while scripts without that extension plots separate points for the samples.

###Analysis of homozygous PTEN deletions, including dilution of tumor purities
Script **PTEN\_CN\_BREASTv7.R** finds the tumor samples in the breast v7 cohort that have a homozygous deletion of tumor suppressor gene PTEN, using WGS data. **plotQDNAseqBREASTv7\_PTENlow\_chr10.R** plots the copy number profiles on chromosome 10 and PTEN for the WGS data for these samples. 

The BAM files with panel data for these samples are investigated in **BAManalysisBv7\_PTENdel.R**. This script performs simulated dilution of the samples with homozygous deletion of PTEN by mixing reads from tumor BAM and normal BAM into a new BAM file. CNVkit is run on panel data for original tumors and normal, and all the dilutions in between, with scripts **runvCNVkitPTENdel\*.R** (early version for original tumor samples: **CNVkitBREASTv7\_PTENdel.txt**). CNVkit with a reference of tumors is run with **plotCompWgsPanelBv7PTENdelTumorRef.R**. Script **workflowPTENdelDilution.R** combines the creation of “diluted” BAM files and run of CNVkit on these in one script. This script is written in a way which only requires change of the desired tumor read content to be adapted for any desired dilution. Reads are counted and added with **bedtoolsCountReads\*.R**/**addReadCounts2Cnr\*.R**/**bedtoolsCountReadsAdd2Cnr\*.R** (however these scripts have been slightly re-written for the breast v6 data). 

Chromosome 10 and PTEN copy number profiles for samples with homozygous PTEN deletion are compared for WGS and panel data with scripts **plotCompWgsPanelBv7PTENdel\*.R**. These scripts also save the data for the bins in PTEN and two control regions, one 5 Mb downstream of PTEN and one 5 Mb upstream of PTEN. p-values between PTEN and control regions are calculated with Mann-Whitney test and saved. The detection capacity of the homozygous deletions of PTEN for different tumor purities is further investigated in script **boxplotBv7PTENdel.R**.

\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

This summary covers most of the scripts in this folder. There are also a few more, which are not described here. These are mostly short and checks some certain feature of the data, but are not used for further analysis. 
