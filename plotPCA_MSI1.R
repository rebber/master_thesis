#Plot PCA components, MAPD vs mrpd and genome-wide segment/bins plots for MSI1 panel data (tumors and normals), treated in CNVkit with either blood reference from 5 normals or a flat reference (log2=0 everywhere)

library(tools)
library(stats)
#install.packages("rgl") 
library(rgl)
library(matrixStats)
library(scatterplot3d)
library(MASS)
path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/"

files_cnr_normals_ref5 <- dir(pattern="txt", paste(path_CNVkit,"batchResultsNormals/cnrRcComb",sep=""), recursive = TRUE)
files_cns_normals_ref5 <- dir(pattern="cns", paste(path_CNVkit,"batchResultsNormals",sep=""), recursive = TRUE)
files_cnr_tumors_ref5 <- dir(pattern="txt", paste(path_CNVkit,"batchResultsRef5/cnrRcComb",sep=""), recursive = TRUE)
#files_cnr_tumors_ref5 <- files_cnr_tumors[-grep("broken_pipe",files_cnr_tumors)]
files_cns_tumors_ref5 <- dir(pattern="cns", paste(path_CNVkit,"batchResultsRef5",sep=""), recursive = TRUE)
#files_cns_tumors_ref5 <- files_cns_tumors[-grep("broken_pipe",files_cns_tumors)]

files_cnr_normals_refFlat <- dir(pattern="B-TD1-CS1-capped.txt", paste(path_CNVkit,"batchResultsRefFlat/cnrRcComb",sep=""), recursive = TRUE)
files_cns_normals_refFlat <- dir(pattern="B-TD1-CS1-capped.cns", paste(path_CNVkit,"batchResultsRefFlat",sep=""), recursive = TRUE)
files_cnr_tumors_refFlat <- dir(pattern="T-TD1-CS1-capped.txt", paste(path_CNVkit,"batchResultsRefFlat/cnrRcComb",sep=""), recursive = TRUE)
#files_cnr_tumors <- files_cnr_tumors[-grep("broken_pipe",files_cnr_tumors)]
files_cns_tumors_refFlat <- dir(pattern="T-TD1-CS1-capped.cns", paste(path_CNVkit,"batchResultsRefFlat",sep=""), recursive = TRUE)
#files_cns_tumors <- files_cns_tumors[-grep("broken_pipe",files_cns_tumors)]

samples <- basename(file_path_sans_ext(c(files_cnr_normals_ref5,files_cnr_tumors_ref5)))

#preallocate error matrix ref 5
example_bins_ref5 <- read.table(paste(path_CNVkit,"batchResultsNormals/cnrRcComb/", files_cnr_normals_ref5[1], sep=""), header = TRUE, as.is=TRUE)
err_mat_ref5 <- matrix(nrow=length(files_cnr_normals_ref5)+length(files_cnr_tumors_ref5),ncol=length(example_bins_ref5$chromosome))
rownames(err_mat_ref5) <- c(files_cnr_normals_ref5,files_cnr_tumors_ref5)

#preallocate error matrix ref flat
example_bins_refFlat <- read.table(paste(path_CNVkit,"batchResultsRefFlat/cnrRcComb/", files_cnr_normals_refFlat[1], sep=""), header = TRUE, as.is=TRUE)
err_mat_refFlat <- matrix(nrow=length(files_cnr_normals_refFlat)+length(files_cnr_tumors_refFlat),ncol=length(example_bins_refFlat$chromosome))
rownames(err_mat_refFlat) <- c(files_cnr_normals_refFlat,files_cnr_tumors_refFlat)

#preallocate for MAPD plot
MAPD_total_ref5 <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(MAPD_total_ref5) <- c("sample","normals","tumors")
MAPD_targets_ref5 <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(MAPD_targets_ref5) <- c("sample","normals","tumors")
MAPD_antitargets_ref5 <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(MAPD_antitargets_ref5) <- c("sample","normals","tumors")
mrpb_total_ref5 <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(mrpb_total_ref5) <- c("sample","normals","tumors")
mrpb_targets_ref5 <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(mrpb_targets_ref5) <- c("sample","normals","tumors")
mrpb_antitargets_ref5 <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(mrpb_antitargets_ref5) <- c("sample","normals","tumors")

MAPD_total_refFlat <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(MAPD_total_refFlat) <- c("sample","normals","tumors")
MAPD_targets_refFlat <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(MAPD_targets_refFlat) <- c("sample","normals","tumors")
MAPD_antitargets_refFlat <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(MAPD_antitargets_refFlat) <- c("sample","normals","tumors")
mrpb_total_refFlat <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(mrpb_total_refFlat) <- c("sample","normals","tumors")
mrpb_targets_refFlat <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(mrpb_targets_refFlat) <- c("sample","normals","tumors")
mrpb_antitargets_refFlat <- matrix(c(samples[1:24],rep(NA,2*length(samples[1:24]))),length(samples[1:24]),3)
colnames(mrpb_antitargets_refFlat) <- c("sample","normals","tumors")

#Get the lengths of all chromosomes and find at which genomic coordinates the chromosomes start and end
chr_sizes <- read.table("/home/rebecka.bergstrom/cnvkit/hg19.chrom.sizes.txt",as.is=TRUE) #length of each chromosome (hg19)
chr_sizes$V1 <- gsub("chr","",chr_sizes$V1)
chr_sizes <- chr_sizes[-grep("_",chr_sizes$V1),]
chr_sizes <- chr_sizes[-grep("M",chr_sizes$V1),]
chr_coord <- as.data.frame(cbind(c(1:22,"X","Y"),rep(NA,24),rep(NA,24)),stringsAsFactors=FALSE) #create matrix with "genome coordinates" for each chromosome (for consecutive chromosomes)
colnames(chr_coord) <-c("chromosome","start","end")
start <- 1 
end <- chr_sizes[which(chr_coord$chromosome[1]==chr_sizes[,1]),2]
chr_coord$start[1] <- start
chr_coord$end[1] <- end
for (i in 2:24) {
  start <- start + chr_sizes[which(chr_coord$chromosome[i-1]==chr_sizes[,1]),2]
  end <- start -1 + chr_sizes[which(chr_coord$chromosome[i]==chr_sizes[,1]),2]
  chr_coord$start[i] <- start
  chr_coord$end[i] <- end
}
chr_coord$start <- as.double(chr_coord$start)

#Find the start and end position of PTEN on chr 10
genes_bed <- read.table("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/QDNAseqWGS/segments2genes/example_file.qdnaseq.genes.bed",as.is=TRUE)
PTEN_ind <- grep("PTEN",genes_bed$V4) #Finds both PTEN and PTENP1
PTEN_ind <- PTEN_ind[-which(PTEN_ind==grep("PTENP1",genes_bed$V4))] 
PTEN_bed <- genes_bed[PTEN_ind,] #The bed-data for PTEN
colnames(PTEN_bed) <- c("chromosome", "start", "end","gene","CN","strand")
PTEN_genestart <- PTEN_bed$start
PTEN_geneend <- PTEN_bed$end
PTEN_genomestart <- as.numeric(chr_coord$end[9]) + PTEN_genestart
PTEN_genomeend <- as.numeric(chr_coord$end[9]) + PTEN_geneend

###########################################
#run for all samples
ptm<-proc.time()
for (i in 1:length(files_cnr_normals_refFlat)) {
  #with reference from 5 MSI1 normals:
  #the data for the normals_ref5
  dat_cnr_normals_ref5 <- read.table(paste(path_CNVkit,"batchResultsNormals/cnrRcComb/", files_cnr_normals_ref5[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_normals_ref5$segLog2 <- rep(NA,length(dat_cnr_normals_ref5$chromosome))
  dat_cnr_normals_ref5$startcoord <- rep(NA,length(dat_cnr_normals_ref5$chromosome))
  dat_cnr_normals_ref5$endcoord <- rep(NA,length(dat_cnr_normals_ref5$chromosome))
  dat_cns_normals_ref5 <- read.table(paste(path_CNVkit,"batchResultsNormals/", files_cns_normals_ref5[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cns_normals_ref5$startcoord <- rep(NA,length(dat_cns_normals_ref5$chromosome))
  dat_cns_normals_ref5$endcoord <- rep(NA,length(dat_cns_normals_ref5$chromosome))
  
  #add the genome wide coordinates
  for (j in 1:24) {
    index_cnr_normals_ref5 <- which(dat_cnr_normals_ref5$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_normals_ref5$startcoord[index_cnr_normals_ref5] <- chr_coord$start[j] + dat_cnr_normals_ref5$start[index_cnr_normals_ref5] - 1 #add the "genome coordinates" for each bin
    dat_cnr_normals_ref5$endcoord[index_cnr_normals_ref5] <- chr_coord$start[j] + dat_cnr_normals_ref5$end[index_cnr_normals_ref5] -1 #add the "genome coordinates" for each bin
    index_cns_normals_ref5 <- which(dat_cns_normals_ref5$chromosome==chr_coord$chromosome[j])
    dat_cns_normals_ref5$startcoord[index_cns_normals_ref5] <- chr_coord$start[j] + dat_cns_normals_ref5$start[index_cns_normals_ref5] - 1 #add the "genome coordinates" for each segment
    dat_cns_normals_ref5$endcoord[index_cns_normals_ref5] <- chr_coord$start[j] + dat_cns_normals_ref5$end[index_cns_normals_ref5] -1 #add the "genome coordinates" for each segment
  }
  
  #add the log2-values of the segment corresponding to each bin
  for (j in 1:length(dat_cns_normals_ref5$chromosome)) {
    ind <- which(dat_cnr_normals_ref5$chromosome==dat_cns_normals_ref5$chromosome[j] & dat_cnr_normals_ref5$start>=dat_cns_normals_ref5$start[j] & dat_cnr_normals_ref5$end<=dat_cns_normals_ref5$end[j])
    dat_cnr_normals_ref5$segLog2[ind] <- dat_cns_normals_ref5$log2[j]
  }
  
  #add the error (log2seg-log2bin) from the normals_ref5 to the error matrix
  err_mat_ref5[i,] <- dat_cnr_normals_ref5$segLog2 - dat_cnr_normals_ref5$log2
  
  #the data for the tumors_ref5
  dat_cnr_tumors_ref5 <- read.table(paste(path_CNVkit,"batchResultsRef5/cnrRcComb/", files_cnr_tumors_ref5[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_tumors_ref5$segLog2 <- rep(NA,length(dat_cnr_tumors_ref5$chromosome))
  dat_cnr_tumors_ref5$startcoord <- rep(NA,length(dat_cnr_tumors_ref5$chromosome))
  dat_cnr_tumors_ref5$endcoord <- rep(NA,length(dat_cnr_tumors_ref5$chromosome))
  dat_cns_tumors_ref5 <- read.table(paste(path_CNVkit,"batchResultsRef5/", files_cns_tumors_ref5[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cns_tumors_ref5$startcoord <- rep(NA,length(dat_cns_tumors_ref5$chromosome))
  dat_cns_tumors_ref5$endcoord <- rep(NA,length(dat_cns_tumors_ref5$chromosome))
  
  #add the genome wide coordinates
  for (j in 1:24) {
    index_cnr_tumors_ref5 <- which(dat_cnr_tumors_ref5$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_tumors_ref5$startcoord[index_cnr_tumors_ref5] <- chr_coord$start[j] + dat_cnr_tumors_ref5$start[index_cnr_tumors_ref5] - 1 #add the "genome coordinates" for each bin
    dat_cnr_tumors_ref5$endcoord[index_cnr_tumors_ref5] <- chr_coord$start[j] + dat_cnr_tumors_ref5$end[index_cnr_tumors_ref5] -1 #add the "genome coordinates" for each bin
    index_cns_tumors_ref5 <- which(dat_cns_tumors_ref5$chromosome==chr_coord$chromosome[j])
    dat_cns_tumors_ref5$startcoord[index_cns_tumors_ref5] <- chr_coord$start[j] + dat_cns_tumors_ref5$start[index_cns_tumors_ref5] - 1 #add the "genome coordinates" for each segment
    dat_cns_tumors_ref5$endcoord[index_cns_tumors_ref5] <- chr_coord$start[j] + dat_cns_tumors_ref5$end[index_cns_tumors_ref5] -1 #add the "genome coordinates" for each segment
  }
  
  #add the log2-values of the segment corresponding to each bin 
  for (j in 1:length(dat_cns_tumors_ref5$chromosome)) {
    ind <- which(dat_cnr_tumors_ref5$chromosome==dat_cns_tumors_ref5$chromosome[j] & dat_cnr_tumors_ref5$start>=dat_cns_tumors_ref5$start[j] & dat_cnr_tumors_ref5$end<=dat_cns_tumors_ref5$end[j])
    dat_cnr_tumors_ref5$segLog2[ind] <- dat_cns_tumors_ref5$log2[j]
  }
  
  #add the error (log2seg-log2bin) from the tumors_ref5 to the error matrix
  err_mat_ref5[length(files_cnr_normals_ref5)+i,] <- dat_cnr_tumors_ref5$segLog2 - dat_cnr_tumors_ref5$log2
  
  #find median read depths per bin (mrpb) and store
  mrpb_total_ref5[i,2] <- median(dat_cnr_normals_ref5$reads)
  mrpb_total_ref5[i,3] <- median(dat_cnr_tumors_ref5$reads)
  mrpb_targets_ref5[i,2] <- median(dat_cnr_normals_ref5$reads[which(dat_cnr_normals_ref5$targetType=="T")])
  mrpb_targets_ref5[i,3] <- median(dat_cnr_tumors_ref5$reads[which(dat_cnr_tumors_ref5$targetType=="T")])
  mrpb_antitargets_ref5[i,2] <- median(dat_cnr_normals_ref5$reads[which(dat_cnr_normals_ref5$targetType=="AT")])
  mrpb_antitargets_ref5[i,3] <- median(dat_cnr_tumors_ref5$reads[which(dat_cnr_tumors_ref5$targetType=="AT")])
  print(paste("The median read depth per bin in sample", samples[i], "for normals_ref5 is", mrpb_targets_ref5[i,2], "reads per bin in targets,", mrpb_antitargets_ref5[i,2], "reads per bin in antitargets and", mrpb_total_ref5[i,2], "reads per bin in total."))
  print(paste("The median read depth per bin in sample", samples[i], "for tumors_ref5 is", mrpb_targets_ref5[i,3], "reads per bin in targets,", mrpb_antitargets_ref5[i,3], "reads per bin in antitargets and",mrpb_total_ref5[i,3], "reads per bin in total."))
  
  #find the median absolute pairwise difference (MAPD) for the log2CN-ratios and store
  MAPD_total_ref5[i,2] <- median(abs(diff(dat_cnr_normals_ref5$log2)))
  MAPD_total_ref5[i,3] <- median(abs(diff(dat_cnr_tumors_ref5$log2)))
  MAPD_targets_ref5[i,2] <- median(abs(diff(dat_cnr_normals_ref5$log2[which(dat_cnr_normals_ref5$targetType=="T")])))
  MAPD_targets_ref5[i,3] <- median(abs(diff(dat_cnr_tumors_ref5$log2[which(dat_cnr_tumors_ref5$targetType=="T")])))
  MAPD_antitargets_ref5[i,2] <- median(abs(diff(dat_cnr_normals_ref5$log2[which(dat_cnr_normals_ref5$targetType=="AT")])))
  MAPD_antitargets_ref5[i,3] <- median(abs(diff(dat_cnr_tumors_ref5$log2[which(dat_cnr_tumors_ref5$targetType=="AT")])))
  
  ###################
  
  #with flat reference:
  #the data for the normals_refFlat
  dat_cnr_normals_refFlat <- read.table(paste(path_CNVkit,"batchResultsRefFlat/cnrRcComb/", files_cnr_normals_refFlat[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_normals_refFlat$segLog2 <- rep(NA,length(dat_cnr_normals_refFlat$chromosome))
  dat_cnr_normals_refFlat$startcoord <- rep(NA,length(dat_cnr_normals_refFlat$chromosome))
  dat_cnr_normals_refFlat$endcoord <- rep(NA,length(dat_cnr_normals_refFlat$chromosome))
  dat_cns_normals_refFlat <- read.table(paste(path_CNVkit,"batchResultsRefFlat/", files_cns_normals_refFlat[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cns_normals_refFlat$startcoord <- rep(NA,length(dat_cns_normals_refFlat$chromosome))
  dat_cns_normals_refFlat$endcoord <- rep(NA,length(dat_cns_normals_refFlat$chromosome))
  
  #add the genome wide coordinates
  for (j in 1:24) {
    index_cnr_normals_refFlat <- which(dat_cnr_normals_refFlat$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_normals_refFlat$startcoord[index_cnr_normals_refFlat] <- chr_coord$start[j] + dat_cnr_normals_refFlat$start[index_cnr_normals_refFlat] - 1 #add the "genome coordinates" for each bin
    dat_cnr_normals_refFlat$endcoord[index_cnr_normals_refFlat] <- chr_coord$start[j] + dat_cnr_normals_refFlat$end[index_cnr_normals_refFlat] -1 #add the "genome coordinates" for each bin
    index_cns_normals_refFlat <- which(dat_cns_normals_refFlat$chromosome==chr_coord$chromosome[j])
    dat_cns_normals_refFlat$startcoord[index_cns_normals_refFlat] <- chr_coord$start[j] + dat_cns_normals_refFlat$start[index_cns_normals_refFlat] - 1 #add the "genome coordinates" for each segment
    dat_cns_normals_refFlat$endcoord[index_cns_normals_refFlat] <- chr_coord$start[j] + dat_cns_normals_refFlat$end[index_cns_normals_refFlat] -1 #add the "genome coordinates" for each segment
  }
  
  #add the log2-values of the segment corresponding to each bin
  for (j in 1:length(dat_cns_normals_refFlat$chromosome)) {
    ind <- which(dat_cnr_normals_refFlat$chromosome==dat_cns_normals_refFlat$chromosome[j] & dat_cnr_normals_refFlat$start>=dat_cns_normals_refFlat$start[j] & dat_cnr_normals_refFlat$end<=dat_cns_normals_refFlat$end[j])
    dat_cnr_normals_refFlat$segLog2[ind] <- dat_cns_normals_refFlat$log2[j]
  }
  
  #add the error (log2seg-log2bin) from the normals_refFlat to the error matrix
  err_mat_refFlat[i,] <- dat_cnr_normals_refFlat$segLog2 - dat_cnr_normals_refFlat$log2
  
  #the data for the tumors_refFlat
  dat_cnr_tumors_refFlat <- read.table(paste(path_CNVkit,"batchResultsRefFlat/cnrRcComb/", files_cnr_tumors_refFlat[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_tumors_refFlat$segLog2 <- rep(NA,length(dat_cnr_tumors_refFlat$chromosome))
  dat_cnr_tumors_refFlat$startcoord <- rep(NA,length(dat_cnr_tumors_refFlat$chromosome))
  dat_cnr_tumors_refFlat$endcoord <- rep(NA,length(dat_cnr_tumors_refFlat$chromosome))
  dat_cns_tumors_refFlat <- read.table(paste(path_CNVkit,"batchResultsRefFlat/", files_cns_tumors_refFlat[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cns_tumors_refFlat$startcoord <- rep(NA,length(dat_cns_tumors_refFlat$chromosome))
  dat_cns_tumors_refFlat$endcoord <- rep(NA,length(dat_cns_tumors_refFlat$chromosome))
  
  #add the genome wide coordinates
  for (j in 1:24) {
    index_cnr_tumors_refFlat <- which(dat_cnr_tumors_refFlat$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_tumors_refFlat$startcoord[index_cnr_tumors_refFlat] <- chr_coord$start[j] + dat_cnr_tumors_refFlat$start[index_cnr_tumors_refFlat] - 1 #add the "genome coordinates" for each bin
    dat_cnr_tumors_refFlat$endcoord[index_cnr_tumors_refFlat] <- chr_coord$start[j] + dat_cnr_tumors_refFlat$end[index_cnr_tumors_refFlat] -1 #add the "genome coordinates" for each bin
    index_cns_tumors_refFlat <- which(dat_cns_tumors_refFlat$chromosome==chr_coord$chromosome[j])
    dat_cns_tumors_refFlat$startcoord[index_cns_tumors_refFlat] <- chr_coord$start[j] + dat_cns_tumors_refFlat$start[index_cns_tumors_refFlat] - 1 #add the "genome coordinates" for each segment
    dat_cns_tumors_refFlat$endcoord[index_cns_tumors_refFlat] <- chr_coord$start[j] + dat_cns_tumors_refFlat$end[index_cns_tumors_refFlat] -1 #add the "genome coordinates" for each segment
  }
  
  #add the log2-values of the segment corresponding to each bin 
  for (j in 1:length(dat_cns_tumors_refFlat$chromosome)) {
    ind <- which(dat_cnr_tumors_refFlat$chromosome==dat_cns_tumors_refFlat$chromosome[j] & dat_cnr_tumors_refFlat$start>=dat_cns_tumors_refFlat$start[j] & dat_cnr_tumors_refFlat$end<=dat_cns_tumors_refFlat$end[j])
    dat_cnr_tumors_refFlat$segLog2[ind] <- dat_cns_tumors_refFlat$log2[j]
  }
  
  #add the error (log2seg-log2bin) from the tumors_refFlat to the error matrix
  err_mat_refFlat[length(files_cnr_normals_refFlat)+i,] <- dat_cnr_tumors_refFlat$segLog2 - dat_cnr_tumors_refFlat$log2
  
  #find median read depths per bin (mrpb) and store
  mrpb_total_refFlat[i,2] <- median(dat_cnr_normals_refFlat$reads)
  mrpb_total_refFlat[i,3] <- median(dat_cnr_tumors_refFlat$reads)
  mrpb_targets_refFlat[i,2] <- median(dat_cnr_normals_refFlat$reads[which(dat_cnr_normals_refFlat$targetType=="T")])
  mrpb_targets_refFlat[i,3] <- median(dat_cnr_tumors_refFlat$reads[which(dat_cnr_tumors_refFlat$targetType=="T")])
  mrpb_antitargets_refFlat[i,2] <- median(dat_cnr_normals_refFlat$reads[which(dat_cnr_normals_refFlat$targetType=="AT")])
  mrpb_antitargets_refFlat[i,3] <- median(dat_cnr_tumors_refFlat$reads[which(dat_cnr_tumors_refFlat$targetType=="AT")])
  print(paste("The median read depth per bin in sample", samples[i], "for normals_refFlat is", mrpb_targets_refFlat[i,2], "reads per bin in targets,", mrpb_antitargets_refFlat[i,2], "reads per bin in antitargets and", mrpb_total_refFlat[i,2], "reads per bin in total."))
  print(paste("The median read depth per bin in sample", samples[i], "for tumors_refFlat is", mrpb_targets_refFlat[i,3], "reads per bin in targets,", mrpb_antitargets_refFlat[i,3], "reads per bin in antitargets and",mrpb_total_refFlat[i,3], "reads per bin in total."))
  
  #find the median absolute pairwise difference (MAPD) for the log2CN-ratios and store
  MAPD_total_refFlat[i,2] <- median(abs(diff(dat_cnr_normals_refFlat$log2)))
  MAPD_total_refFlat[i,3] <- median(abs(diff(dat_cnr_tumors_refFlat$log2)))
  MAPD_targets_refFlat[i,2] <- median(abs(diff(dat_cnr_normals_refFlat$log2[which(dat_cnr_normals_refFlat$targetType=="T")])))
  MAPD_targets_refFlat[i,3] <- median(abs(diff(dat_cnr_tumors_refFlat$log2[which(dat_cnr_tumors_refFlat$targetType=="T")])))
  MAPD_antitargets_refFlat[i,2] <- median(abs(diff(dat_cnr_normals_refFlat$log2[which(dat_cnr_normals_refFlat$targetType=="AT")])))
  MAPD_antitargets_refFlat[i,3] <- median(abs(diff(dat_cnr_tumors_refFlat$log2[which(dat_cnr_tumors_refFlat$targetType=="AT")])))
  
  #################
  #plot genomewide and PTEN
  setwd(paste(path_CNVkit,"Rplots/refComp/",sep=""))
  
  #the indeces and coordinates for zoom in on PTEN 
  startZoomIn <- PTEN_genomestart - 50*15000 #genome coordinate of the first bin included in zoom in
  endZoomIn <- PTEN_genomeend + 50*15000 #genome coordinate of last bin included in zoom in
  #for the normals with ref 5
  zoomIn_cnr_normals_ref5 <- which(dat_cnr_normals_ref5$startcoord>startZoomIn & dat_cnr_normals_ref5$startcoord<endZoomIn) #which bins start in the zoom in window
  zoomIn_cns_normals_ref5 <- which(dat_cns_normals_ref5$startcoord>startZoomIn & dat_cns_normals_ref5$startcoord<endZoomIn) #which segments start in the zoom in window
  segmStartXcoordZoomIn_normals_ref5 <- c(startZoomIn,dat_cns_normals_ref5$startcoord[zoomIn_cns_normals_ref5]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn_normals_ref5 <- c(dat_cns_normals_ref5$endcoord[zoomIn_cns_normals_ref5-1],endZoomIn) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns_normals_ref5)>0) {
    segmYcoordZoomIn_normals_ref5 <- dat_cns_normals_ref5$log2[c((zoomIn_cns_normals_ref5[1]-1),zoomIn_cns_normals_ref5)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn_normals_ref5 <- dat_cns_normals_ref5$log2[which(dat_cns_normals_ref5$startcoord<=startZoomIn & dat_cns_normals_ref5$endcoord>=startZoomIn)]
  }
  #for the tumors with ref5
  zoomIn_cnr_tumors_ref5 <- which(dat_cnr_tumors_ref5$startcoord>startZoomIn & dat_cnr_tumors_ref5$startcoord<endZoomIn) #which bins start in the zoom in window
  zoomIn_cns_tumors_ref5 <- which(dat_cns_tumors_ref5$startcoord>startZoomIn & dat_cns_tumors_ref5$startcoord<endZoomIn) #which segments start in the zoom in window
  segmStartXcoordZoomIn_tumors_ref5 <- c(startZoomIn,dat_cns_tumors_ref5$startcoord[zoomIn_cns_tumors_ref5]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn_tumors_ref5 <- c(dat_cns_tumors_ref5$endcoord[zoomIn_cns_tumors_ref5-1],endZoomIn) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns_tumors_ref5)>0) {
    segmYcoordZoomIn_tumors_ref5 <- dat_cns_tumors_ref5$log2[c((zoomIn_cns_tumors_ref5[1]-1),zoomIn_cns_tumors_ref5)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn_tumors_ref5 <- dat_cns_tumors_ref5$log2[which(dat_cns_tumors_ref5$startcoord<=startZoomIn & dat_cns_tumors_ref5$endcoord>=startZoomIn)]
  }
  #for the normals with flat ref
  zoomIn_cnr_normals_refFlat <- which(dat_cnr_normals_refFlat$startcoord>startZoomIn & dat_cnr_normals_refFlat$startcoord<endZoomIn) #which bins start in the zoom in window
  zoomIn_cns_normals_refFlat <- which(dat_cns_normals_refFlat$startcoord>startZoomIn & dat_cns_normals_refFlat$startcoord<endZoomIn) #which segments start in the zoom in window
  segmStartXcoordZoomIn_normals_refFlat <- c(startZoomIn,dat_cns_normals_refFlat$startcoord[zoomIn_cns_normals_refFlat]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn_normals_refFlat <- c(dat_cns_normals_refFlat$endcoord[zoomIn_cns_normals_refFlat-1],endZoomIn) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns_normals_refFlat)>0) {
    segmYcoordZoomIn_normals_refFlat <- dat_cns_normals_refFlat$log2[c((zoomIn_cns_normals_refFlat[1]-1),zoomIn_cns_normals_refFlat)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn_normals_refFlat <- dat_cns_normals_refFlat$log2[which(dat_cns_normals_refFlat$startcoord<=startZoomIn & dat_cns_normals_refFlat$endcoord>=startZoomIn)]
  }
  #for the tumors with flat ref
  zoomIn_cnr_tumors_refFlat <- which(dat_cnr_tumors_refFlat$startcoord>startZoomIn & dat_cnr_tumors_refFlat$startcoord<endZoomIn) #which bins start in the zoom in window
  zoomIn_cns_tumors_refFlat <- which(dat_cns_tumors_refFlat$startcoord>startZoomIn & dat_cns_tumors_refFlat$startcoord<endZoomIn) #which segments start in the zoom in window
  segmStartXcoordZoomIn_tumors_refFlat <- c(startZoomIn,dat_cns_tumors_refFlat$startcoord[zoomIn_cns_tumors_refFlat]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn_tumors_refFlat <- c(dat_cns_tumors_refFlat$endcoord[zoomIn_cns_tumors_refFlat-1],endZoomIn) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns_tumors_refFlat)>0) {
    segmYcoordZoomIn_tumors_refFlat <- dat_cns_tumors_refFlat$log2[c((zoomIn_cns_tumors_refFlat[1]-1),zoomIn_cns_tumors_refFlat)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn_tumors_refFlat <- dat_cns_tumors_refFlat$log2[which(dat_cns_tumors_refFlat$startcoord<=startZoomIn & dat_cns_tumors_refFlat$endcoord>=startZoomIn)]
  }
  
  #plot genome wide comparison between wgs data and panel data with tumor reference
  title1 <- paste("log2 copy number CNVkit, ref 5 blood, MSI1 normal, ",basename(file_path_sans_ext(files_cnr_normals_ref5[i])),sep="")
  title2 <- paste("log2 copy number CNVkit, ref flat, MSI1 normal, ",basename(file_path_sans_ext(files_cnr_normals_refFlat[i])),sep="")
  title3 <- paste("log2 copy number CNVkit, ref 5 blood, MSI1 tumor, ",basename(file_path_sans_ext(files_cnr_tumors_ref5[i])),sep="")
  title4 <- paste("log2 copy number CNVkit, ref flat, MSI1 tumor, ",basename(file_path_sans_ext(files_cnr_tumors_refFlat[i])),sep="")
  title_zoomIn <- "zoom in on PTEN"
  
  filename1 <- paste("MSI1_compRefs_genome_PTEN_",basename(file_path_sans_ext(gsub("-TD1-CS1-capped","",files_cnr_normals_ref5[i]))),".pdf",sep="")
  filename2 <- paste("MSI1_compRefs_genome_PTEN_",basename(file_path_sans_ext(gsub("-TD1-CS1-capped","",files_cnr_tumors_ref5[i]))),".pdf",sep="")
  xticks <- as.numeric(chr_coord$start)+(as.numeric(chr_coord$end)-as.numeric(chr_coord$start))/2
  
  #plot for normals with ref 5 and ref flat in the same graph
  pdf(filename1,width=20,height=9.8)
  layout(matrix(c(1,1,2,3,3,4), 2, 3, byrow = TRUE)) 
  #plot normals ref 5 
  plot(dat_cnr_normals_ref5$startcoord, dat_cnr_normals_ref5$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title1,xlab="Chromosome #", ylab="log2 CN ratio")
  segments(dat_cns_normals_ref5$startcoord, dat_cns_normals_ref5$log2, dat_cns_normals_ref5$endcoord, dat_cns_normals_ref5$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"))
  abline(v=PTEN_genomestart,lty=2)
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  #zoom in on PTEN 
  plot(dat_cnr_normals_ref5$startcoord[zoomIn_cnr_normals_ref5], dat_cnr_normals_ref5$log2[zoomIn_cnr_normals_ref5], xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN ratio")
  segments(segmStartXcoordZoomIn_normals_ref5,segmYcoordZoomIn_normals_ref5,segmEndXcoordZoomIn_normals_ref5,segmYcoordZoomIn_normals_ref5, lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7))
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  #plot normals ref flat 
  plot(dat_cnr_normals_refFlat$startcoord, dat_cnr_normals_refFlat$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title2,xlab="Chromosome #", ylab="log2 CN ratio")
  segments(dat_cns_normals_refFlat$startcoord, dat_cns_normals_refFlat$log2, dat_cns_normals_refFlat$endcoord, dat_cns_normals_refFlat$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"))
  abline(v=PTEN_genomestart,lty=2)
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  #zoom in on PTEN 
  plot(dat_cnr_normals_refFlat$startcoord[zoomIn_cnr_normals_refFlat], dat_cnr_normals_refFlat$log2[zoomIn_cnr_normals_refFlat], xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN ratio")
  segments(segmStartXcoordZoomIn_normals_refFlat,segmYcoordZoomIn_normals_refFlat,segmEndXcoordZoomIn_normals_refFlat,segmYcoordZoomIn_normals_refFlat, lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7))
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  dev.off()
  
  #plot for tumors with ref 5 and ref flat in the same graph
  pdf(filename2,width=20,height=9.8)
  layout(matrix(c(1,1,2,3,3,4), 2, 3, byrow = TRUE)) 
  #plot tumors ref 5 
  plot(dat_cnr_tumors_ref5$startcoord, dat_cnr_tumors_ref5$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title3,xlab="Chromosome #", ylab="log2 CN ratio")
  segments(dat_cns_tumors_ref5$startcoord, dat_cns_tumors_ref5$log2, dat_cns_tumors_ref5$endcoord, dat_cns_tumors_ref5$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"))
  abline(v=PTEN_genomestart,lty=2)
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  #zoom in on PTEN 
  plot(dat_cnr_tumors_ref5$startcoord[zoomIn_cnr_tumors_ref5], dat_cnr_tumors_ref5$log2[zoomIn_cnr_tumors_ref5], xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN ratio")
  segments(segmStartXcoordZoomIn_tumors_ref5,segmYcoordZoomIn_tumors_ref5,segmEndXcoordZoomIn_tumors_ref5,segmYcoordZoomIn_tumors_ref5, lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7))
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  #plot tumors ref flat 
  plot(dat_cnr_tumors_refFlat$startcoord, dat_cnr_tumors_refFlat$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title4,xlab="Chromosome #", ylab="log2 CN ratio")
  segments(dat_cns_tumors_refFlat$startcoord, dat_cns_tumors_refFlat$log2, dat_cns_tumors_refFlat$endcoord, dat_cns_tumors_refFlat$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"))
  abline(v=PTEN_genomestart,lty=2)
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  #zoom in on PTEN 
  plot(dat_cnr_tumors_refFlat$startcoord[zoomIn_cnr_tumors_refFlat], dat_cnr_tumors_refFlat$log2[zoomIn_cnr_tumors_refFlat], xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN ratio")
  segments(segmStartXcoordZoomIn_tumors_refFlat,segmYcoordZoomIn_tumors_refFlat,segmEndXcoordZoomIn_tumors_refFlat,segmYcoordZoomIn_tumors_refFlat, lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7))
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  dev.off()
  
}
proc.time()-ptm

#count bins in PTEN:
length(which(dat_cnr_tumors_refFlat$startcoord>PTEN_genomestart & dat_cnr_tumors_refFlat$startcoord<PTEN_genomeend & dat_cnr_tumors_refFlat$targetType=="T"))

#tables of coverage
View(cbind(mrpb_total_ref5,mrpb_targets_ref5[,2:3],mrpb_antitargets_ref5[,2:3],signif(as.numeric(mrpb_targets_ref5[,2])/as.numeric(mrpb_antitargets_ref5[,2]),digits=3),signif(as.numeric(mrpb_targets_ref5[,3])/as.numeric(mrpb_antitargets_ref5[,3]),digits=3)))
colMedians(cbind(mrpb_total_ref5,mrpb_targets_ref5[,2:3],mrpb_antitargets_ref5[,2:3],signif(as.numeric(mrpb_targets_ref5[,2])/as.numeric(mrpb_antitargets_ref5[,2]),digits=3),signif(as.numeric(mrpb_targets_ref5[,3])/as.numeric(mrpb_antitargets_ref5[,3]),digits=3)))
View(cbind(mrpb_total_refFlat,mrpb_targets_refFlat[,2:3],mrpb_antitargets_refFlat[,2:3],signif(as.numeric(mrpb_targets_refFlat[,2])/as.numeric(mrpb_antitargets_refFlat[,2]),digits=3),signif(as.numeric(mrpb_targets_refFlat[,3])/as.numeric(mrpb_antitargets_refFlat[,3]),digits=3)))

#######################################

#plot MAPD plots
setwd(paste(path_CNVkit,"Rplots/",sep=""))
reads=1:2000
mapd<-c()
for (i in 1:2000) mapd[i]=median(abs(diff(log2(rpois(10000,reads[i])/reads[i]))))
filename_MAPD <- "MAPD_MSI1_NormTum.pdf"
pdf(filename_MAPD,width=9,height=9)
layout(matrix(1,1,1,byrow=TRUE))
plot(reads,mapd,xlab='median #reads per bin',ylab='median absolute pairwise difference (log2CNratio)',type='l',main='Minimum stochastic noise, MSI1',log="x",ylim=c(0,1))
points(reads, mapd*sqrt(2),type="l",lty=2)
#points(as.double(mrpb_total_ref5[,2]),MAPD_total_ref5[,2],col="red")
#points(as.double(mrpb_total_ref5[,3]),MAPD_total_ref5[,3],col="blue")
cols <- c(topo.colors(8))
points(as.double(mrpb_targets_ref5[,2]),MAPD_targets_ref5[,2],col=cols[1],pch=16)
points(as.double(mrpb_targets_ref5[,3]),MAPD_targets_ref5[,3],col=cols[2],pch=16)
points(as.double(mrpb_antitargets_ref5[,2]),MAPD_antitargets_ref5[,2],col=cols[3],pch=16)
points(as.double(mrpb_antitargets_ref5[,3]),MAPD_antitargets_ref5[,3],col=cols[4],pch=16)
points(as.double(mrpb_targets_refFlat[,2]),MAPD_targets_refFlat[,2],col=cols[5],pch=16)
points(as.double(mrpb_targets_refFlat[,3]),MAPD_targets_refFlat[,3],col=cols[6],pch=16)
points(as.double(mrpb_antitargets_refFlat[,2]),MAPD_antitargets_refFlat[,2],col=cols[7],pch=16)
points(as.double(mrpb_antitargets_refFlat[,3]),MAPD_antitargets_refFlat[,3],col=cols[8],pch=16)
legend("topright", c("targets, normals, ref5","targets, tumors, ref5","antitargets, normals, ref5","antitargets, tumors, ref5","targets, normals, refFlat","targets, tumors, refFlat","antitargets, normals, refFlat","antitargets, tumors, refFlat","Poisson distribution","Poisson*sqrt(2)"),col=c(cols,"black","black"),pch=c(rep(16,8),NA,NA),lty=c(rep(NA,8),1,2))
dev.off()

layout(matrix(1,1,1,byrow=TRUE))
plot(reads,mapd,xlab='median #reads per bin',ylab='median absolute pairwise difference (log2CNratio)',type='l',main='Minimum stochastic noise, MSI1',log="x",ylim=c(0,1))
points(reads, mapd*sqrt(2),type="l",lty=2)
#points(as.double(mrpb_total_ref5[,2]),MAPD_total_ref5[,2],col="red")
#points(as.double(mrpb_total_ref5[,3]),MAPD_total_ref5[,3],col="blue")
cols <- "black" #c(topo.colors(8))
text(x=as.double(mrpb_targets_ref5[,2]),y=MAPD_targets_ref5[,2],col=cols[1],labels=as.character(1:24))
text(as.double(mrpb_targets_ref5[,3]),MAPD_targets_ref5[,3],col=cols[2],labels=as.character(25:48))
text(as.double(mrpb_antitargets_ref5[,2]),MAPD_antitargets_ref5[,2],col=cols[3],labels=as.character(1:24))
text(as.double(mrpb_antitargets_ref5[,3]),MAPD_antitargets_ref5[,3],col=cols[4],labels=as.character(25:48))
text(as.double(mrpb_targets_refFlat[,2]),MAPD_targets_refFlat[,2],col=cols[5],labels=as.character(1:24))
text(as.double(mrpb_targets_refFlat[,3]),MAPD_targets_refFlat[,3],col=cols[6],labels=as.character(25:48))
text(as.double(mrpb_antitargets_refFlat[,2]),MAPD_antitargets_refFlat[,2],col=cols[7],labels=as.character(1:24))
text(as.double(mrpb_antitargets_refFlat[,3]),MAPD_antitargets_refFlat[,3],col=cols[8],labels=as.character(25:48))
legend("topright", c("targets, normals, ref5","targets, tumors, ref5","antitargets, normals, ref5","antitargets, tumors, ref5","targets, normals, refFlat","targets, tumors, refFlat","antitargets, normals, refFlat","antitargets, tumors, refFlat","Poisson distribution","Poisson*sqrt(2)"),col=c(cols,"black","black"),pch=c(rep(16,8),NA,NA),lty=c(rep(NA,8),1,2))


#########################

#PCA

#ref of 5 normals
#since some bins are not covered by a segment their error is NA, which s not tolerated by princomp()
#therefore exchange the NAs to 0 (the errors should be centered around 0 and just a few extra bins with this value should not be a problem)
#(another alternative would be to replace the NAs with the respective median value for that certain bin, but I don't do that now)
ind_na_ref5 <- which(is.na(err_mat_ref5))
err_mat_ref5[ind_na_ref5] <- 0
#err_mat_ref5[ind_na_ref5] <- NA

pca_ref5 <- prcomp(err_mat_ref5)
#pca_ref5 <- prcomp(t(err_mat_ref5))
#pca_ref5 <- princomp(t(err_mat_ref5))
plot(pca_ref5)
plot(pca_ref5$x)
box_ref5 <- boxplot(pca_ref5$x,plot=FALSE)
pca_ref5$x[,1]
out1_ref5 <- box_ref5$out[which(box_ref5$group==1)]
for (i in 1:length(out1_ref5)) {
  print(paste("Sample",samples[which(pca_ref5$x[,1]==out1_ref5[i])],"is outlier with value:",signif(pca_ref5$x[which(pca_ref5$x[,1]==out1_ref5[i]),1],digits=3)))
}
quantile(pca_ref5$x[,1])
pca_ref5$x[,2]
out2_ref5 <- box_ref5$out[which(box_ref5$group==2)]
for (i in 1:length(out2_ref5)) {
  print(paste("Sample",samples[which(pca_ref5$x[,2]==out2_ref5[i])],"is outlier with value:",signif(pca_ref5$x[which(pca_ref5$x[,2]==out2_ref5[i]),1],digits=3)))
}
quantile(pca_ref5$x[,2])
plot(pca_ref5$x[,3],pca_ref5$x[,4])
plot(pca_ref5$x[,4],pca_ref5$x[,5])

#remove samples with failed segmentation from pca analysis
rm_samples <- c(grep("rem1T",samples),grep("rem2T",samples),grep("rem3T",samples))
#exchange "rem1" etc for correct sample IDs
err_mat_ref5_red <- err_mat_ref5[-rm_samples,]
pca_ref5_red <- prcomp(err_mat_ref5_red)
plot(pca_ref5_red$x[,1],pca_ref5_red$x[,2])
plot(pca_ref5_red$x[25:45,1],pca_ref5_red$x[25:45,2],col="red")
points(pca_ref5_red$x[1:24,1],pca_ref5_red$x[1:24,2])
plot(pca_ref5_red$x[,1],pca_ref5_red$x[,2],type="n",ylim=c(-35,35),xlim=c(-50,0))
text(pca_ref5_red$x[,1],pca_ref5_red$x[,2],pch=as.character(1:45))
plot(pca_ref5_red$x[,2],pca_ref5_red$x[,3])
plot(pca_ref5_red$x[,3],pca_ref5_red$x[,4])

#remove three further samples with mrpb=<100
rm_samples2 <- c(rm_samples,grep("rem4T",samples),grep("rem5T",samples),grep("rem6T",samples))
#exchange "rem1" etc for correct sample IDs
err_mat_ref5_red2 <- err_mat_ref5[-rm_samples2,]
pca_ref5_red2 <- prcomp(err_mat_ref5_red2)
layout(matrix(1,1,1))
plot(pca_ref5_red2$x[,1],pca_ref5_red2$x[,2])
plot(pca_ref5_red2$x[25:42,1],pca_ref5_red2$x[25:42,2],col="red",ylim=c(-300,300))
points(pca_ref5_red2$x[1:24,1],pca_ref5_red2$x[1:24,2])
plot(pca_ref5_red2$x[,1],pca_ref5_red2$x[,2],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_ref5_red2$x[,1],pca_ref5_red2$x[,2],pch=as.character(1:42))
plot(pca_ref5_red2$x[,2],pca_ref5_red2$x[,3],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_ref5_red2$x[,2],pca_ref5_red2$x[,3],pch=as.character(1:42))
plot(pca_ref5_red2$x[,3],pca_ref5_red2$x[,4],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_ref5_red2$x[,3],pca_ref5_red2$x[,4],pch=as.character(1:42))
plot(pca_ref5_red2$x[,2],pca_ref5_red2$x[,3])
plot(pca_ref5_red2$x[,3],pca_ref5_red2$x[,4])

#remove another three samples with mrpb<200 (in total 9 samples removed)
rm_samples3 <- c(rm_samples2,grep("rem7T",samples),grep("rem8T",samples),grep("rem9T",samples))
#exchange "rem1" etc for correct sample IDs
err_mat_ref5_red3 <- err_mat_ref5[-rm_samples3,]
pca_ref5_red3 <- prcomp(err_mat_ref5_red3)
layout(matrix(1,1,1))
plot(pca_ref5_red3$x[,1],pca_ref5_red3$x[,2])
plot(pca_ref5_red3$x[25:38,1],pca_ref5_red3$x[25:38,2],col="red",ylim=c(-300,300))
points(pca_ref5_red3$x[1:24,1],pca_ref5_red3$x[1:24,2])
plot(pca_ref5_red3$x[,1],pca_ref5_red3$x[,2],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_ref5_red3$x[,1],pca_ref5_red3$x[,2],pch=as.character(1:42))
plot(pca_ref5_red3$x[,2],pca_ref5_red3$x[,3],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_ref5_red3$x[,2],pca_ref5_red3$x[,3],pch=as.character(1:42))
plot(pca_ref5_red3$x[,3],pca_ref5_red3$x[,4],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_ref5_red3$x[,3],pca_ref5_red3$x[,4],pch=as.character(1:42))
plot(pca_ref5_red3$x[,2],pca_ref5_red3$x[,3])
plot(pca_ref5_red3$x[,3],pca_ref5_red3$x[,4])

#flat ref
ind_na_refFlat <- which(is.na(err_mat_refFlat))
err_mat_refFlat[ind_na_refFlat] <- 0
#err_mat_refFlat[ind_na_refFlat] <- NA

pca_refFlat <- prcomp(err_mat_refFlat)
#pca_refFlat <- prcomp(t(err_mat_refFlat))
#pca_refFlat <- princomp(t(err_mat_refFlat))
plot(pca_refFlat)
plot(pca_refFlat$x)
box_refFlat <- boxplot(pca_refFlat$x,plot=FALSE)
pca_refFlat$x[,1]
out1_refFlat <- box_refFlat$out[which(box_refFlat$group==1)]
for (i in 1:length(out1_refFlat)) {
  print(paste("Sample",samples[which(pca_refFlat$x[,1]==out1_refFlat[i])],"is outlier with value:",signif(pca_refFlat$x[which(pca_refFlat$x[,1]==out1_refFlat[i]),1],digits=3)))
}
quantile(pca_refFlat$x[,1])
pca_refFlat$x[,2]
out2_refFlat <- box_refFlat$out[which(box_refFlat$group==2)]
for (i in 1:length(out2_refFlat)) {
  print(paste("Sample",samples[which(pca_refFlat$x[,2]==out2_refFlat[i])],"is outlier with value:",signif(pca_refFlat$x[which(pca_refFlat$x[,2]==out2_refFlat[i]),1],digits=3)))
}
quantile(pca_refFlat$x[,2])
plot(pca_refFlat$x[,3],pca_refFlat$x[,4])
plot(pca_refFlat$x[,4],pca_refFlat$x[,5])

#remove samples with failed segmentation from pca analysis
rm_samples <- c(grep("rem1T",samples),grep("rem2T",samples),grep("rem3T",samples))
#exchange "rem1" etc for correct sample IDs
err_mat_refFlat_red <- err_mat_refFlat[-rm_samples,]
pca_refFlat_red <- prcomp(err_mat_refFlat_red)
plot(pca_refFlat_red$x[,1],pca_refFlat_red$x[,2])
plot(pca_refFlat_red$x[25:45,1],pca_refFlat_red$x[25:45,2],col="red")
points(pca_refFlat_red$x[1:24,1],pca_refFlat_red$x[1:24,2])
plot(pca_refFlat_red$x[,1],pca_refFlat_red$x[,2],type="n",ylim=c(-35,35),xlim=c(-50,0))
text(pca_refFlat_red$x[,1],pca_refFlat_red$x[,2],pch=as.character(1:45))
plot(pca_refFlat_red$x[,2],pca_refFlat_red$x[,3])
plot(pca_refFlat_red$x[,3],pca_refFlat_red$x[,4])
View(samples[-rm_samples])

#remove three further samples with mrpb=<100
rm_samples2 <- c(rm_samples,grep("rem4T",samples),grep("rem5T",samples),grep("rem6T",samples))
#exchange "rem1" etc for correct sample IDs
err_mat_refFlat_red2 <- err_mat_refFlat[-rm_samples2,]
pca_refFlat_red2 <- prcomp(err_mat_refFlat_red2)
layout(matrix(1,1,1))
plot(pca_refFlat_red2$x[,1],pca_refFlat_red2$x[,2])
plot(pca_refFlat_red2$x[25:42,1],pca_refFlat_red2$x[25:42,2],col="red",ylim=c(-300,300))
points(pca_refFlat_red2$x[1:24,1],pca_refFlat_red2$x[1:24,2])
plot(pca_refFlat_red2$x[,1],pca_refFlat_red2$x[,2],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_refFlat_red2$x[,1],pca_refFlat_red2$x[,2],pch=as.character(1:42))
plot(pca_refFlat_red2$x[,2],pca_refFlat_red2$x[,3],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_refFlat_red2$x[,2],pca_refFlat_red2$x[,3],pch=as.character(1:42))
plot(pca_refFlat_red2$x[,3],pca_refFlat_red2$x[,4],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_refFlat_red2$x[,3],pca_refFlat_red2$x[,4],pch=as.character(1:42))
plot(pca_refFlat_red2$x[,2],pca_refFlat_red2$x[,3])
plot(pca_refFlat_red2$x[,3],pca_refFlat_red2$x[,4])

#remove another three samples with mrpb<200 (in total 9 samples removed)
rm_samples3 <- c(rm_samples2,grep("rem7T",samples),grep("rem8T",samples),grep("rem9T",samples))
#exchange "rem1" etc for correct sample IDs
err_mat_refFlat_red3 <- err_mat_refFlat[-rm_samples3,]
pca_refFlat_red3 <- prcomp(err_mat_refFlat_red3)
layout(matrix(1,1,1))
plot(pca_refFlat_red3$x[,1],pca_refFlat_red3$x[,2])
plot(pca_refFlat_red3$x[25:38,1],pca_refFlat_red3$x[25:38,2],col="red",ylim=c(-300,300))
points(pca_refFlat_red3$x[1:24,1],pca_refFlat_red3$x[1:24,2])
plot(pca_refFlat_red3$x[,1],pca_refFlat_red3$x[,2],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_refFlat_red3$x[,1],pca_refFlat_red3$x[,2],pch=as.character(1:42))
plot(pca_refFlat_red3$x[,2],pca_refFlat_red3$x[,3],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_refFlat_red3$x[,2],pca_refFlat_red3$x[,3],pch=as.character(1:42))
plot(pca_refFlat_red3$x[,3],pca_refFlat_red3$x[,4],type="n") #,ylim=c(-35,35),xlim=c(-50,0)
text(pca_refFlat_red3$x[,3],pca_refFlat_red3$x[,4],pch=as.character(1:42))
plot(pca_refFlat_red3$x[,2],pca_refFlat_red3$x[,3])
plot(pca_refFlat_red3$x[,3],pca_refFlat_red3$x[,4])

#save pca-plots
setwd(paste(path_CNVkit,"Rplots/pca",sep=""))
filename_pca <- "pca_MSI1_allSamples.pdf"
pdf(filename_pca,width=10,height=9)
layout(matrix(1,1,1,byrow=TRUE))
plot(pca_ref5$x[,1],pca_ref5$x[,2],xlab="PC1",ylab="PC2",main="PCA, all MSI1 samples, 5 blood samples as ref")
plot(pca_ref5$x[,3],pca_ref5$x[,4],xlab="PC3",ylab="PC4",main="PCA, all MSI1 samples, 5 blood samples as ref")
plot(pca_ref5$x[,4],pca_ref5$x[,5],xlab="PC4",ylab="PC5",main="PCA, all MSI1 samples, 5 blood samples as ref")
plot(pca_refFlat$x[,1],pca_refFlat$x[,2],xlab="PC1",ylab="PC2",main="PCA, all MSI1 samples, flat ref")
plot(pca_refFlat$x[,3],pca_refFlat$x[,4],xlab="PC3",ylab="PC4",main="PCA, all MSI1 samples, flat ref")
plot(pca_refFlat$x[,4],pca_refFlat$x[,5],xlab="PC4",ylab="PC5",main="PCA, all MSI1 samples, flat ref")
dev.off()

filename_pca_red <- "pca_MSI1_no3lowCov.pdf"
pdf(filename_pca_red,width=15,height=9)
layout(matrix(1,1,1,byrow=TRUE))
plot(pca_ref5_red$x[,1],pca_ref5_red$x[,2],xlab="PC1",ylab="PC2",main="PCA, MSI1 samples (not 3 low cov), 5 blood samples as ref")
plot(pca_ref5_red$x[25:45,1],pca_ref5_red$x[25:45,2],col="red",xlab="PC1",ylab="PC2",main="PCA, MSI1 samples (not 3 low cov), 5 blood samples as ref")
points(pca_ref5_red$x[1:24,1],pca_ref5_red$x[1:24,2],col="blue")
legend("topleft",c("Normals","Tumors"),pch=c(1,1),col=c("blue","red"))
plot(pca_ref5_red$x[,2],pca_ref5_red$x[,3],xlab="PC2",ylab="PC3",main="PCA, MSI1 samples (not 3 low cov), 5 blood samples as ref")
plot(pca_ref5_red$x[,3],pca_ref5_red$x[,4],xlab="PC3",ylab="PC4",main="PCA, MSI1 samples (not 3 low cov), 5 blood samples as ref")
plot(pca_refFlat_red$x[,1],pca_refFlat_red$x[,2],xlab="PC1",ylab="PC2",main="PCA, MSI1 samples (not 3 low cov), flat ref")
plot(pca_refFlat_red$x[25:45,1],pca_refFlat_red$x[25:45,2],col="red",xlab="PC1",ylab="PC2",main="PCA, MSI1 samples (not 3 low cov), flat ref")
points(pca_refFlat_red$x[1:24,1],pca_refFlat_red$x[1:24,2],col="blue")
legend("topleft",c("Normals","Tumors"),pch=c(1,1),col=c("blue","red"))
plot(pca_refFlat_red$x[,2],pca_refFlat_red$x[,3],xlab="PC2",ylab="PC3",main="PCA, MSI1 samples (not 3 low cov), flat ref")
plot(pca_refFlat_red$x[,3],pca_refFlat_red$x[,4],xlab="PC3",ylab="PC4",main="PCA, MSI1 samples (not 3 low cov), flat ref")
dev.off()


#############

#MDS (same as PCA)
mds_ref5_red <- cmdscale(dist(err_mat_ref5_red),k=44)
plot(mds_ref5_red[,1],mds_ref5_red[,2])
plot(mds_ref5_red[,2],mds_ref5_red[,3])
plot(mds_ref5_red[,3],mds_ref5_red[,4])

#isoMDS (gives also the same as PCA)
imds_ref5_red <- isoMDS(dist(err_mat_ref5_red),k=44)
plot(imds_ref5_red$points[,1],imds_ref5_red$points[,2])
plot(imds_ref5_red$points[,2],imds_ref5_red$points[,3])
plot(imds_ref5_red$points[,3],imds_ref5_red$points[,4])

###########

#plot error for each target bin towards pc 1:8 when 3 samples are removed
targ <- which(dat_cnr_tumors_ref5$targetType=="T")
err_mat_ref5_red_targ <- err_mat_ref5_red[,targ]
for(i in 1:length(err_mat_ref5_red[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca_ref5_red$x[,1],err_mat_ref5_red_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red$x[,1],err_mat_ref5_red_targ[,i],label=1:45)
  plot(pca_ref5_red$x[,2],err_mat_ref5_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_ref5_red$x[,2],err_mat_ref5_red_targ[,i],label=1:45)
  plot(pca_ref5_red$x[,3],err_mat_ref5_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_ref5_red$x[,3],err_mat_ref5_red_targ[,i],label=1:45)
  plot(pca_ref5_red$x[,4],err_mat_ref5_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_ref5_red$x[,4],err_mat_ref5_red_targ[,i],label=1:45)
  plot(pca_ref5_red$x[,5],err_mat_ref5_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_ref5_red$x[,5],err_mat_ref5_red_targ[,i],label=1:45)
  plot(pca_ref5_red$x[,6],err_mat_ref5_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_ref5_red$x[,6],err_mat_ref5_red_targ[,i],label=1:45)
  plot(pca_ref5_red$x[,7],err_mat_ref5_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_ref5_red$x[,7],err_mat_ref5_red_targ[,i],label=1:45)
  plot(pca_ref5_red$x[,8],err_mat_ref5_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_ref5_red$x[,8],err_mat_ref5_red_targ[,i],label=1:45)
  browser()
}

#when 6 samples are removed
err_mat_ref5_red2_targ <- err_mat_ref5_red2[,targ]
for(i in 1:length(err_mat_ref5_red2[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca_ref5_red2$x[,1],err_mat_ref5_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red2$x[,1],err_mat_ref5_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_ref5_red2$x[,2],err_mat_ref5_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red2$x[,2],err_mat_ref5_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_ref5_red2$x[,3],err_mat_ref5_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red2$x[,3],err_mat_ref5_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_ref5_red2$x[,4],err_mat_ref5_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red2$x[,4],err_mat_ref5_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_ref5_red2$x[,5],err_mat_ref5_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red2$x[,5],err_mat_ref5_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_ref5_red2$x[,6],err_mat_ref5_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red2$x[,6],err_mat_ref5_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_ref5_red2$x[,7],err_mat_ref5_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red2$x[,7],err_mat_ref5_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_ref5_red2$x[,8],err_mat_ref5_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red2$x[,8],err_mat_ref5_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  browser()
}

#when 9 samples are removed
err_mat_ref5_red3_targ <- err_mat_ref5_red3[,targ]
for(i in 1:length(err_mat_ref5_red3[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca_ref5_red3$x[,1],err_mat_ref5_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red3$x[,1],err_mat_ref5_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_ref5_red3$x[,2],err_mat_ref5_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red3$x[,2],err_mat_ref5_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_ref5_red3$x[,3],err_mat_ref5_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red3$x[,3],err_mat_ref5_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_ref5_red3$x[,4],err_mat_ref5_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red3$x[,4],err_mat_ref5_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_ref5_red3$x[,5],err_mat_ref5_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red3$x[,5],err_mat_ref5_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_ref5_red3$x[,6],err_mat_ref5_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red3$x[,6],err_mat_ref5_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_ref5_red3$x[,7],err_mat_ref5_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red3$x[,7],err_mat_ref5_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_ref5_red3$x[,8],err_mat_ref5_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_ref5_red3$x[,8],err_mat_ref5_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  browser()
}

#for refFlat:
#plot error for each target bin towards pc 1:8 when 3 samples are removed
targ <- which(dat_cnr_tumors_refFlat$targetType=="T")
err_mat_refFlat_red_targ <- err_mat_refFlat_red[,targ]
for(i in 1:length(err_mat_refFlat_red[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca_refFlat_red$x[,1],err_mat_refFlat_red_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red$x[,1],err_mat_refFlat_red_targ[,i],label=1:45)
  plot(pca_refFlat_red$x[,2],err_mat_refFlat_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_refFlat_red$x[,2],err_mat_refFlat_red_targ[,i],label=1:45)
  plot(pca_refFlat_red$x[,3],err_mat_refFlat_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_refFlat_red$x[,3],err_mat_refFlat_red_targ[,i],label=1:45)
  plot(pca_refFlat_red$x[,4],err_mat_refFlat_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_refFlat_red$x[,4],err_mat_refFlat_red_targ[,i],label=1:45)
  plot(pca_refFlat_red$x[,5],err_mat_refFlat_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_refFlat_red$x[,5],err_mat_refFlat_red_targ[,i],label=1:45)
  plot(pca_refFlat_red$x[,6],err_mat_refFlat_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_refFlat_red$x[,6],err_mat_refFlat_red_targ[,i],label=1:45)
  plot(pca_refFlat_red$x[,7],err_mat_refFlat_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_refFlat_red$x[,7],err_mat_refFlat_red_targ[,i],label=1:45)
  plot(pca_refFlat_red$x[,8],err_mat_refFlat_red_targ[,i],ylim=c(-1.5,1.5),type="n")
  text(pca_refFlat_red$x[,8],err_mat_refFlat_red_targ[,i],label=1:45)
  browser()
}

#when 6 samples are removed
err_mat_refFlat_red2_targ <- err_mat_refFlat_red2[,targ]
for(i in 1:length(err_mat_refFlat_red2[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca_refFlat_red2$x[,1],err_mat_refFlat_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red2$x[,1],err_mat_refFlat_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_refFlat_red2$x[,2],err_mat_refFlat_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red2$x[,2],err_mat_refFlat_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_refFlat_red2$x[,3],err_mat_refFlat_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red2$x[,3],err_mat_refFlat_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_refFlat_red2$x[,4],err_mat_refFlat_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red2$x[,4],err_mat_refFlat_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_refFlat_red2$x[,5],err_mat_refFlat_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red2$x[,5],err_mat_refFlat_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_refFlat_red2$x[,6],err_mat_refFlat_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red2$x[,6],err_mat_refFlat_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_refFlat_red2$x[,7],err_mat_refFlat_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red2$x[,7],err_mat_refFlat_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  plot(pca_refFlat_red2$x[,8],err_mat_refFlat_red2_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red2$x[,8],err_mat_refFlat_red2_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples2]))
  browser()
}

#when 9 samples are removed
err_mat_refFlat_red3_targ <- err_mat_refFlat_red3[,targ]
for(i in 1:length(err_mat_refFlat_red3[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca_refFlat_red3$x[,1],err_mat_refFlat_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red3$x[,1],err_mat_refFlat_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_refFlat_red3$x[,2],err_mat_refFlat_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red3$x[,2],err_mat_refFlat_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_refFlat_red3$x[,3],err_mat_refFlat_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red3$x[,3],err_mat_refFlat_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_refFlat_red3$x[,4],err_mat_refFlat_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red3$x[,4],err_mat_refFlat_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_refFlat_red3$x[,5],err_mat_refFlat_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red3$x[,5],err_mat_refFlat_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_refFlat_red3$x[,6],err_mat_refFlat_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red3$x[,6],err_mat_refFlat_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_refFlat_red3$x[,7],err_mat_refFlat_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red3$x[,7],err_mat_refFlat_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  plot(pca_refFlat_red3$x[,8],err_mat_refFlat_red3_targ[,i],ylim=c(-1.5,1.5),type="n",main=paste("i=",i))
  text(pca_refFlat_red3$x[,8],err_mat_refFlat_red3_targ[,i],label=gsub("-TD1-CS1-capped","",samples[-rm_samples3]))
  browser()
}

############################

#plot 3d-plot with MAPD (targets), mrpb (targets) and pc1 for tumors (run code in terminal to get window with movable 3d plot)
plot3d(mrpb_targets_ref5[-c(rm_samples-24),3],MAPD_targets_ref5[-c(rm_samples-24),3],pca_ref5_red$x[25:45,1],main="Ref 5")
plot3d(mrpb_targets_refFlat[-c(rm_samples-24),3],MAPD_targets_refFlat[-c(rm_samples-24),3],pca_refFlat_red$x[25:25,1],main="Ref flat")
plot3d(mrpb_targets_ref5[-c(rm_samples2-24),3],MAPD_targets_ref5[-c(rm_samples2-24),3],pca_ref5_red2$x[25:42,1],main="Ref 5")
plot3d(mrpb_targets_ref5[-c(rm_samples3 -24),3],MAPD_targets_ref5[-c(rm_samples3-24),3],pca_ref5_red3$x[25:39,1],main="Ref 5")
plot3d(mrpb_targets_refFlat[-c(rm_samples-24),3],MAPD_targets_refFlat[-c(rm_samples-24),3],pca_refFlat_red$x[25:45,1],main="Ref flat")
plot3d(mrpb_targets_refFlat[-c(rm_samples2-24),3],MAPD_targets_refFlat[-c(rm_samples2-24),3],pca_refFlat_red2$x[25:42,1],main="Ref flat")
plot3d(mrpb_targets_refFlat[-c(rm_samples3-24),3],MAPD_targets_refFlat[-c(rm_samples3-24),3],pca_refFlat_red3$x[25:39,1],main="Ref flat")

setwd(paste(path_CNVkit,"Rplots",sep=""))
pdf("3d_MAPD_mrpb_pc1.pdf")
layout(matrix(1,1,1,byrow=TRUE))
scatterplot3d(mrpb_targets_ref5[-c(rm_samples-24,25:48),3],MAPD_targets_ref5[-c(rm_samples-24,25:48),3],pca_ref5_red$x[25:45,1],color=topo.colors(21),pch=16,highlight.3d = TRUE)
dev.off()

plot(pca_ref5_red$x[25:45,1],MAPD_targets_ref5[-c(rm_samples-24,25:48),3])
plot(pca_ref5_red$x[25:45,6],MAPD_targets_ref5[-c(rm_samples-24,25:48),3])

plot(mrpb_targets_ref5[-c(25:48),3],pca_ref5$x[25:48,1])
plot(mrpb_targets_ref5[-c(rm_samples-24,25:48),3],pca_ref5_red$x[25:45,1])





