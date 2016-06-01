#Plot PCA components, MAPD and mrpd for Bv6 panel data (tumors), treated in CNVkit with a flat reference (log2=0 everywhere)
#Created by Rebecka Bergström on 2016-03-24

#load the PCA workspace created with this script
#load("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/Rplots/smooth/randTestSet/pca/withoutY/.RData")


library(tools)
library(stats)
#install.packages("rgl") 
library(rgl)
library(matrixStats)
library(scatterplot3d)
library(MASS)
library(permute)
library(dbscan)
library(class)

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/"

#find the files
files_cnr_tumors_refFlat <- dir(pattern="T*panel_v1.txt", paste(path_CNVkit,"cnrRcComb/addedSmooth/",sep=""), recursive = TRUE)
files_cns_tumors_refFlat <- dir(pattern="T*panel_v1.cns", paste(path_CNVkit,"batchResults",sep=""), recursive = TRUE)
files_cns_tumors_refFlat <- append(files_cns_tumors_refFlat,dir(pattern="T*panel_v1.cns", paste(path_CNVkit,"batchResultsDevs",sep=""), recursive = TRUE))
files_cns_tumors_refFlat <- files_cns_tumors_refFlat[order(files_cns_tumors_refFlat)]
samples <- basename(file_path_sans_ext(files_cnr_tumors_refFlat))

#load variable testset, defining which samples belong to the test set
load(paste(path_CNVkit,"test50/test50",sep=""))

#preallocate error matrix ref flat
example_bins_refFlat <- read.table(paste(path_CNVkit,"cnrRcComb/addedSmooth/", files_cnr_tumors_refFlat[1], sep=""), header = TRUE, as.is=TRUE)
err_mat_refFlat <- matrix(nrow=length(samples),ncol=length(example_bins_refFlat$chromosome))
rownames(err_mat_refFlat) <- c(files_cnr_tumors_refFlat)

#preallocate matrix with smoothed values for each bin
smooth_mat_refFlat <- matrix(nrow=length(samples),ncol=length(example_bins_refFlat$chromosome))
rownames(smooth_mat_refFlat) <- c(files_cnr_tumors_refFlat)

#preallocate for MAPD plot
MAPD_total_refFlat <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(MAPD_total_refFlat) <- c("sample","tumors tot")
MAPD_targets_refFlat <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(MAPD_targets_refFlat) <- c("sample","tumors targ")
MAPD_antitargets_refFlat <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),)
colnames(MAPD_antitargets_refFlat) <- c("sample","tumors atarg")
mrpb_total_refFlat <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(mrpb_total_refFlat) <- c("sample","tumors tot")
mrpb_targets_refFlat <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(mrpb_targets_refFlat) <- c("sample","tumors targ")
mrpb_antitargets_refFlat <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(mrpb_antitargets_refFlat) <- c("sample","tumors atarg")

#preallocate tot_reads matrix
tot_reads <- matrix(rep(NA,2*length(samples)),length(samples),2)
colnames(tot_reads) <- c("total_rc","onTargetRate")
rownames(tot_reads) <- samples

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
for (i in 1:length(files_cnr_tumors_refFlat)) { 
  ptm2<-proc.time()
  #with flat reference:
  #the data for the tumors_refFlat
  dat_cnr_tumors_refFlat <- read.table(paste(path_CNVkit,"cnrRcComb/addedSmooth/", files_cnr_tumors_refFlat[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_tumors_refFlat$startcoord <- rep(NA,length(dat_cnr_tumors_refFlat$chromosome))
  dat_cnr_tumors_refFlat$endcoord <- rep(NA,length(dat_cnr_tumors_refFlat$chromosome))
   
  #add the genome wide coordinates
  for (j in 1:24) {
    index_cnr_tumors_refFlat <- which(dat_cnr_tumors_refFlat$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_tumors_refFlat$startcoord[index_cnr_tumors_refFlat] <- chr_coord$start[j] + dat_cnr_tumors_refFlat$start[index_cnr_tumors_refFlat] - 1 #add the "genome coordinates" for each bin
    dat_cnr_tumors_refFlat$endcoord[index_cnr_tumors_refFlat] <- chr_coord$start[j] + dat_cnr_tumors_refFlat$end[index_cnr_tumors_refFlat] -1 #add the "genome coordinates" for each bin
  }
  
  #add the total read count and target read part to the tot_reads matrix
  tot_reads[i,1] <- sum(dat_cnr_tumors_refFlat$reads)
  tot_reads[i,2] <- sum(dat_cnr_tumors_refFlat$reads[which(dat_cnr_tumors_refFlat$targetType=="T")])/tot_reads[i,1]
  
  #add the error (log2bin-smooth) from the tumors_refFlat to the error matrix
  err_mat_refFlat[i,] <- dat_cnr_tumors_refFlat$log2 - dat_cnr_tumors_refFlat$smooth
  
  #add the smooth values for each bin to the smooth matrix
  smooth_mat_refFlat[i,] <- dat_cnr_tumors_refFlat$smooth
  
  #find median read depths per bin (mrpb) and store
  mrpb_total_refFlat[i,2] <- median(dat_cnr_tumors_refFlat$reads)
  mrpb_targets_refFlat[i,2] <- median(dat_cnr_tumors_refFlat$reads[which(dat_cnr_tumors_refFlat$targetType=="T")])
  mrpb_antitargets_refFlat[i,2] <- median(dat_cnr_tumors_refFlat$reads[which(dat_cnr_tumors_refFlat$targetType=="AT")])
  print(paste("The median read depth per bin in sample", samples[i], "for tumors_refFlat is", mrpb_targets_refFlat[i,2], "reads per bin in targets,", mrpb_antitargets_refFlat[i,2], "reads per bin in antitargets and",mrpb_total_refFlat[i,2], "reads per bin in total."))
  
  #find the median absolute pairwise difference (MAPD) for the log2CN-ratios and store
  MAPD_total_refFlat[i,2] <- median(abs(diff(dat_cnr_tumors_refFlat$log2)))
  MAPD_targets_refFlat[i,2] <- median(abs(diff(dat_cnr_tumors_refFlat$log2[which(dat_cnr_tumors_refFlat$targetType=="T")])))
  MAPD_antitargets_refFlat[i,2] <- median(abs(diff(dat_cnr_tumors_refFlat$log2[which(dat_cnr_tumors_refFlat$targetType=="AT")])))
  
  #################
  #plot genomewide and PTEN
  setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/smoothsBinsGenomeWide/",sep=""))
  
  #the indeces and coordinates for zoom in on PTEN 
  startZoomIn <- PTEN_genomestart - 50*15000 #genome coordinate of the first bin included in zoom in
  endZoomIn <- PTEN_genomeend + 50*15000 #genome coordinate of last bin included in zoom in
  #for the tumors with flat ref
  zoomIn_cnr_tumors_refFlat <- which(dat_cnr_tumors_refFlat$startcoord>startZoomIn & dat_cnr_tumors_refFlat$startcoord<endZoomIn) #which bins start in the zoom in window
  
  #plot genome wide comparison between wgs data and panel data with tumor reference
  title4 <- paste("log2 copy number CNVkit, ref flat, Bv6 tumor, ",basename(file_path_sans_ext(files_cnr_tumors_refFlat[i])),sep="")
  title_zoomIn <- "zoom in on PTEN"
  filename2 <- paste("Bv6_genome_PTEN_",basename(file_path_sans_ext(gsub("-TD1-CS1-capped","",files_cnr_tumors_refFlat[i]))),".jpeg",sep="")
  xticks <- as.numeric(chr_coord$start)+(as.numeric(chr_coord$end)-as.numeric(chr_coord$start))/2
  
  #plot for tumors with ref flat 
  jpeg(filename2,width=2000,height=900)
  layout(matrix(c(1,1,2), 1, 3, byrow = TRUE)) 
  #plot tumors ref flat 
  plot(dat_cnr_tumors_refFlat$startcoord, dat_cnr_tumors_refFlat$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title4,xlab="Chromosome #", ylab="log2 CN ratio")
  segments(dat_cnr_tumors_refFlat$startcoord, dat_cnr_tumors_refFlat$smooth, dat_cnr_tumors_refFlat$endcoord, dat_cnr_tumors_refFlat$smooth, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"))
  abline(v=PTEN_genomestart,lty=2)
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  #zoom in on PTEN 
  plot(dat_cnr_tumors_refFlat$startcoord[zoomIn_cnr_tumors_refFlat], dat_cnr_tumors_refFlat$log2[zoomIn_cnr_tumors_refFlat], xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN ratio")
  segments(dat_cnr_tumors_refFlat$startcoord[zoomIn_cnr_tumors_refFlat], dat_cnr_tumors_refFlat$smooth[zoomIn_cnr_tumors_refFlat], dat_cnr_tumors_refFlat$endcoord[zoomIn_cnr_tumors_refFlat], dat_cnr_tumors_refFlat$smooth[zoomIn_cnr_tumors_refFlat], lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7))
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  dev.off()
  proc.time()-ptm2
}
proc.time()-ptm
#end of for-loop for the 334 tumors


#count bins in PTEN:
length(which(dat_cnr_tumors_refFlat$startcoord>PTEN_genomestart & dat_cnr_tumors_refFlat$startcoord<PTEN_genomeend & dat_cnr_tumors_refFlat$targetType=="T"))

#tables of coverage
mrpb_summary <- cbind(mrpb_total_refFlat,mrpb_targets_refFlat[,2],mrpb_antitargets_refFlat[,2],signif(as.numeric(mrpb_targets_refFlat[,2])/as.numeric(mrpb_antitargets_refFlat[,2]),digits=3))
colnames(mrpb_summary) <- c("sample","tumors_tot", "targets", "antitargets", "ratio_t/at")
mrpb_summary <- rbind(mrpb_summary, c("medians",as.character(apply(matrix(as.double(mrpb_summary[,2:5]),length(mrpb_summary[,1]),4),2,median))))
write.table(mrpb_summary,file=paste(path_CNVkit,"Rplots/smooth/randTestSet/MAPD/mrpb_pca_refFlat.txt",sep=""),quote=FALSE,row.names=FALSE)

MAPD_summary <- cbind(MAPD_total_refFlat,MAPD_targets_refFlat[,2],MAPD_antitargets_refFlat[,2],signif(as.numeric(MAPD_targets_refFlat[,2])/as.numeric(MAPD_antitargets_refFlat[,2]),digits=3))
colnames(MAPD_summary) <- c("sample","tumors_tot", "targets", "antitargets", "ratio_t/at")
MAPD_summary <- rbind(MAPD_summary, c("medians",as.character(apply(matrix(as.double(MAPD_summary[,2:5]),length(MAPD_summary[,1]),4),2,median))))
write.table(MAPD_summary,file=paste(path_CNVkit,"Rplots/smooth/randTestSet/MAPD/MAPD_pca_refFlat.txt",sep=""),quote=FALSE,row.names=FALSE)

#######################################

#plot MAPD plots
setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/MAPD",sep=""))
reads=1:7500
mapd<-c()
for (i in 1:7500) mapd[i]=median(abs(diff(log2(rpois(10000,reads[i])/reads[i]))))
filename_MAPD <- "MAPD_Bv6_Tum.pdf"
pdf(filename_MAPD,width=9,height=9)
layout(matrix(1,1,1,byrow=TRUE))
plot(reads,mapd,xlab='median #reads per bin',ylab='median absolute pairwise difference (log2CNratio)',type='l',main='Minimum stochastic noise, Bv6 all 334 tumors',log="x",ylim=c(0,1))
points(reads, mapd*sqrt(2),type="l",lty=2)
#points(as.double(mrpb_total_ref5[,2]),MAPD_total_ref5[,2],col="red")
#points(as.double(mrpb_total_ref5[,3]),MAPD_total_ref5[,3],col="blue")
cols <- c("red","blue")
points(as.double(mrpb_targets_refFlat[,2]),MAPD_targets_refFlat[,2],col=cols[1],pch=16)
points(as.double(mrpb_antitargets_refFlat[,2]),MAPD_antitargets_refFlat[,2],col=cols[2],pch=16)
legend("topright", c("targets, tumors, flat ref","antitargets, tumors, flat ref","Poisson distribution","Poisson*sqrt(2)"),col=c(cols,"black","black"),pch=c(rep(16,2),NA,NA),lty=c(rep(NA,2),1,2))
dev.off()

#########################

#PCA

#flat ref
#remove all bins in Y chromosome since they confuse the PCA
#err_mat_refFlat_total <- err_mat_refFlat #backup
Ybins <- which(dat_cnr_tumors_refFlat$startcoord >= chr_coord$start[24])
#err_mat_refFlat <- err_mat_refFlat_total[,-Ybins]

#do PCA on the training set (i.e. all samples except the 50 in test set)
pca_refFlat <- prcomp(err_mat_refFlat[-as.integer(testset[,2]),-Ybins])
layout(matrix(1,1,1))
barplot(pca_refFlat$sdev)

layout(matrix(c(1:4),2,2, byrow = TRUE))
par(mar=c(5,4,4,2)+0.1)
plot(pca_refFlat$x[,1],pca_refFlat$x[,2])
plot(pca_refFlat$x[,2],pca_refFlat$x[,3])
plot(pca_refFlat$x[,3],pca_refFlat$x[,4])
plot(pca_refFlat$x[,4],pca_refFlat$x[,5])

layout(matrix(c(1:4),2,2, byrow = TRUE))
plot(pca_refFlat$x[,5],pca_refFlat$x[,6])
plot(pca_refFlat$x[,6],pca_refFlat$x[,7])
plot(pca_refFlat$x[,7],pca_refFlat$x[,8])
plot(pca_refFlat$x[,8],pca_refFlat$x[,9])


#save pca-plots
setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/",sep=""))
filename_pca <- "pca_Bv6_allTumors_woY.pdf"
pdf(filename_pca,width=9.8,height=9.8)
layout(matrix(c(1:4),2,2, byrow = TRUE))
plot(pca_refFlat$x[,1],pca_refFlat$x[,2],xlab="PC1",ylab="PC2",main="PCA, all Bv6 tumors, flat ref",xlim=c(-130,130),ylim=c(-130,130))
plot(pca_refFlat$x[,2],pca_refFlat$x[,3],xlab="PC2",ylab="PC3",main="PCA, all Bv6 tumors, flat ref",xlim=c(-130,130),ylim=c(-130,130))
plot(pca_refFlat$x[,3],pca_refFlat$x[,4],xlab="PC3",ylab="PC4",main="PCA, all Bv6 tumors, flat ref",xlim=c(-130,130),ylim=c(-130,130))
plot(pca_refFlat$x[,4],pca_refFlat$x[,5],xlab="PC4",ylab="PC5",main="PCA, all Bv6 tumors, flat ref",xlim=c(-130,130),ylim=c(-130,130))
layout(matrix(c(1:4),2,2, byrow = TRUE))
plot(pca_refFlat$x[,5],pca_refFlat$x[,6],xlab="PC5",ylab="PC6",main="PCA, all Bv6 tumors, flat ref",xlim=c(-130,130),ylim=c(-130,130))
plot(pca_refFlat$x[,6],pca_refFlat$x[,7],xlab="PC6",ylab="PC7",main="PCA, all Bv6 tumors, flat ref",xlim=c(-130,130),ylim=c(-130,130))
plot(pca_refFlat$x[,7],pca_refFlat$x[,8],xlab="PC7",ylab="PC8",main="PCA, all Bv6 tumors, flat ref",xlim=c(-130,130),ylim=c(-130,130))
plot(pca_refFlat$x[,8],pca_refFlat$x[,9],xlab="PC8",ylab="PC9",main="PCA, all Bv6 tumors, flat ref",xlim=c(-130,130),ylim=c(-130,130))
dev.off()


###########

#plot error for each target bin towards pc 1:8 
targ <- which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="T") 
err_mat_refFlat_targ <- err_mat_refFlat[,targ]
for(i in 1:length(err_mat_refFlat[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca_refFlat$x[,1],err_mat_refFlat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n",main=paste("target bin ",i))
  text(pca_refFlat$x[,1],err_mat_refFlat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca_refFlat$x[,2],err_mat_refFlat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,2],err_mat_refFlat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca_refFlat$x[,3],err_mat_refFlat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,3],err_mat_refFlat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca_refFlat$x[,4],err_mat_refFlat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,4],err_mat_refFlat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca_refFlat$x[,5],err_mat_refFlat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,5],err_mat_refFlat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca_refFlat$x[,6],err_mat_refFlat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,6],err_mat_refFlat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca_refFlat$x[,7],err_mat_refFlat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,7],err_mat_refFlat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca_refFlat$x[,8],err_mat_refFlat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,8],err_mat_refFlat_targ[-as.integer(testset[,2]),i],label=1:281)
  browser()
}


######################

#PCA for "negative control"
#shuffle samples within each bin -> same stdev as before
err_mat_shuff <- err_mat_refFlat[-as.integer(testset[,2]),]
ptm <- proc.time()
err_mat_shuff <- apply(err_mat_shuff,2,sample)
proc.time()-ptm
pca_shuff <- prcomp(err_mat_shuff[,-Ybins])
layout(matrix(1,1,1))
barplot(pca_shuff$sdev)
layout(matrix(c(1:4),2,2, byrow = TRUE))
plot(pca_shuff$x[,1],pca_shuff$x[,2],xlab="PC1",ylab="PC2",main="PCA, all Bv6 tumors, flat ref")
plot(pca_shuff$x[,2],pca_shuff$x[,3],xlab="PC2",ylab="PC3",main="PCA, all Bv6 tumors, flat ref")
plot(pca_shuff$x[,3],pca_shuff$x[,4],xlab="PC3",ylab="PC4",main="PCA, all Bv6 tumors, flat ref")
plot(pca_shuff$x[,4],pca_shuff$x[,5],xlab="PC4",ylab="PC5",main="PCA, all Bv6 tumors, flat ref")
layout(matrix(c(1:4),2,2, byrow = TRUE))
plot(pca_shuff$x[,5],pca_shuff$x[,6],xlab="PC5",ylab="PC6",main="PCA, all Bv6 tumors, flat ref")
plot(pca_shuff$x[,6],pca_shuff$x[,7],xlab="PC6",ylab="PC7",main="PCA, all Bv6 tumors, flat ref")
plot(pca_shuff$x[,7],pca_shuff$x[,8],xlab="PC7",ylab="PC8",main="PCA, all Bv6 tumors, flat ref")
plot(pca_shuff$x[,8],pca_shuff$x[,9],xlab="PC8",ylab="PC9",main="PCA, all Bv6 tumors, flat ref")

targ <- which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="T")
err_mat_shuff_targ <- err_mat_shuff[,targ]
for(i in 1:length(err_mat_shuff[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca_shuff$x[,1],err_mat_shuff_targ[,i],ylim=c(-1,1),type="n",main=paste("target bin ",i))
  text(pca_shuff$x[,1],err_mat_shuff_targ[,i],label=1:281)
  plot(pca_shuff$x[,2],err_mat_shuff_targ[,i],ylim=c(-1,1),type="n")
  text(pca_shuff$x[,2],err_mat_shuff_targ[,i],label=1:281)
  plot(pca_shuff$x[,3],err_mat_shuff_targ[,i],ylim=c(-1,1),type="n")
  text(pca_shuff$x[,3],err_mat_shuff_targ[,i],label=1:281)
  plot(pca_shuff$x[,4],err_mat_shuff_targ[,i],ylim=c(-1,1),type="n")
  text(pca_shuff$x[,4],err_mat_shuff_targ[,i],label=1:281)
  plot(pca_shuff$x[,5],err_mat_shuff_targ[,i],ylim=c(-1,1),type="n")
  text(pca_shuff$x[,5],err_mat_shuff_targ[,i],label=1:281)
  plot(pca_shuff$x[,6],err_mat_shuff_targ[,i],ylim=c(-1,1),type="n")
  text(pca_shuff$x[,6],err_mat_shuff_targ[,i],label=1:281)
  plot(pca_shuff$x[,7],err_mat_shuff_targ[,i],ylim=c(-1,1),type="n")
  text(pca_shuff$x[,7],err_mat_shuff_targ[,i],label=1:281)
  plot(pca_shuff$x[,8],err_mat_shuff_targ[,i],ylim=c(-1,1),type="n")
  text(pca_shuff$x[,8],err_mat_shuff_targ[,i],label=1:281)
  browser()
}

#compare our data and control
layout(1,1,1)
barplot(pca_refFlat$sdev-pca_shuff$sdev)
plot(pca_shuff$sdev,pca_refFlat$sdev)

#save plots
setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/",sep=""))
filename_negctrl <- "pca_compNegCtrl_Bv6_woY.pdf"
pdf(filename_negctrl, width=12, height=9.8)
layout(matrix(c(1,2),1,2,byrow=TRUE))
par(mar=c(5,4,4,2)+0.1)
barplot(pca_refFlat$sdev,ylim=c(0,30),xlab="PC",ylab="Standard deviation",main="SD per PC, pca on Bv6 tumors")
barplot(pca_shuff$sdev,ylim=c(0,30),xlab="PC",ylab="Standard deviation",main="SD per PC, pca on neg ctrl")
layout(1,1,1)
plot(pca_shuff$sdev,pca_refFlat$sdev,xlab="SD per PC, neg ctrl",ylab="SD per PC, Bv6 tumors",main="SD per PC, Bv6 tumors vs neg ctrl")
fit<-rlm(pca_refFlat$sdev ~ pca_shuff$sdev)
points(pca_shuff$sdev,fit$fitted.values,type="l")
legend("topleft",c("PCs",paste("Fitted line, slope = ",signif(fit$coefficients[[2]],digits=3),", interc = ",signif(fit$coefficients[[1]],digits=3),sep="")),lty=c(NA,1),pch=c(1,NA),col=1)
dev.off()


######################

#Compare with smoothed values
#PCA on smooth mat
#remove Y bins
smooth_mat_refFlat_total <- smooth_mat_refFlat
#smooth_mat_refFlat <- smooth_mat_refFlat_total[,-Ybins]

pca_smooth <- prcomp(smooth_mat_refFlat[-as.integer(testset[,2]),-Ybins])
layout(matrix(1:8,2,4,byrow=TRUE))
plot(pca_refFlat$x[,1],pca_refFlat$x[,2],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_refFlat$x[,1],pca_refFlat$x[,2])
plot(pca_refFlat$x[,2],pca_refFlat$x[,3],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_refFlat$x[,2],pca_refFlat$x[,3])
plot(pca_refFlat$x[,3],pca_refFlat$x[,4],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_refFlat$x[,3],pca_refFlat$x[,4])
plot(pca_refFlat$x[,4],pca_refFlat$x[,5],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_refFlat$x[,4],pca_refFlat$x[,5])
plot(pca_smooth$x[,1],pca_smooth$x[,2],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_smooth$x[,1],pca_smooth$x[,2])
plot(pca_smooth$x[,2],pca_smooth$x[,3],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_smooth$x[,2],pca_smooth$x[,3])
plot(pca_smooth$x[,3],pca_smooth$x[,4],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_smooth$x[,3],pca_smooth$x[,4])
plot(pca_smooth$x[,4],pca_smooth$x[,5],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_smooth$x[,4],pca_smooth$x[,5])

cor_err_smooth <-cor(pca_refFlat$x,pca_smooth$x) #columns refer to PCs in pca_smooth$x, while rows refer to PCs in pca_refFlat$x
quantile(abs(cor_err_smooth),probs=seq(0,1,0.1))
quantile(abs(cor_err_smooth[1:5,]),probs=seq(0,1,0.1)) #look at PC 1-5 in pca_refFlat$x (the interesting PCs) 
which(cor_err_smooth==max(cor_err_smooth[1:5,]))
length(which(abs(cor_err_smooth)>0.1))
length(which(abs(cor_err_smooth[1:5,])>0.1))
length(which(abs(cor_err_smooth[1,])>0.1))
length(which(abs(cor_err_smooth[2,])>0.1))
which(abs(cor_err_smooth[1,])>0.1)
which(abs(cor_err_smooth[2,])>0.1)
which(abs(cor_err_smooth[3,])>0.1)
which(abs(cor_err_smooth[4,])>0.1)
which(abs(cor_err_smooth[5,])>0.1)
layout(1,1,1)
plot(which(abs(cor_err_smooth[1,])>0.1),rep(1,length(which(abs(cor_err_smooth[1,])>0.1))),ylim=c(0,6),xlim=c(0,400),ylab="PC for errors",xlab="PCs for smooth with cor>0.1")
points(which(abs(cor_err_smooth[2,])>0.1),rep(2,length(which(abs(cor_err_smooth[2,])>0.1))))
points(which(abs(cor_err_smooth[3,])>0.1),rep(3,length(which(abs(cor_err_smooth[3,])>0.1))))
points(which(abs(cor_err_smooth[4,])>0.1),rep(4,length(which(abs(cor_err_smooth[4,])>0.1))))
points(which(abs(cor_err_smooth[5,])>0.1),rep(5,length(which(abs(cor_err_smooth[5,])>0.1))))

length(which(abs(cor_err_smooth)>0.2))
length(which(abs(cor_err_smooth[1:5,])>0.2))
which(abs(cor_err_smooth[1,])>0.2)
which(abs(cor_err_smooth[2,])>0.2)
which(abs(cor_err_smooth[3,])>0.2)
which(abs(cor_err_smooth[4,])>0.2)
which(abs(cor_err_smooth[5,])>0.2)
plot(which(abs(cor_err_smooth[1,])>0.2),rep(1,length(which(abs(cor_err_smooth[1,])>0.2))),ylim=c(0,6),xlim=c(0,50),ylab="PC for errors",xlab="PCs for smooth with cor>0.2")
points(which(abs(cor_err_smooth[2,])>0.2),rep(2,length(which(abs(cor_err_smooth[2,])>0.2))))
points(which(abs(cor_err_smooth[3,])>0.2),rep(3,length(which(abs(cor_err_smooth[3,])>0.2))))
points(which(abs(cor_err_smooth[4,])>0.2),rep(4,length(which(abs(cor_err_smooth[4,])>0.2))))
points(which(abs(cor_err_smooth[5,])>0.2),rep(5,length(which(abs(cor_err_smooth[5,])>0.2))))

layout(matrix(1:10,2,5,byrow=TRUE))
plot(pca_smooth$x[,1],pca_refFlat$x[,1])
plot(pca_smooth$x[,1],pca_refFlat$x[,2])
plot(pca_smooth$x[,1],pca_refFlat$x[,3])
plot(pca_smooth$x[,1],pca_refFlat$x[,4])
plot(pca_smooth$x[,1],pca_refFlat$x[,5])
plot(pca_smooth$x[,2],pca_refFlat$x[,1])
plot(pca_smooth$x[,2],pca_refFlat$x[,2])
plot(pca_smooth$x[,2],pca_refFlat$x[,3])
plot(pca_smooth$x[,2],pca_refFlat$x[,4])
plot(pca_smooth$x[,2],pca_refFlat$x[,5])

which(abs(cor_err_smooth[1:5,])>0.2)
cor_err_smooth[1:5,]

median(cor_err_smooth[1:4,6])
median(cor_err_smooth[1:4,7])
median(cor_err_smooth[1:4,8])
median(cor_err_smooth[1:4,9])
median(cor_err_smooth[1:4,10])


#Shuffle segments values randomly
smooth_mat_shuff <- smooth_mat_refFlat[-as.integer(testset[,2]),]
smooth_mat_shuff <- apply(smooth_mat_shuff,2,sample)
pca_smooth_shuff <- prcomp(smooth_mat_shuff[,-Ybins])
plot(pca_smooth_shuff$sdev,pca_smooth$sdev)
barplot(pca_smooth$sdev-pca_smooth_shuff$sdev)
barplot(pca_smooth$sdev)
barplot(pca_smooth_shuff$sdev)


#compare error-pca and smooth-pca including Ybins
pca_refFlat_total <- prcomp(err_mat_refFlat[-as.integer(testset[,2]),])
pca_smooth_total <- prcomp(smooth_mat_refFlat[-as.integer(testset[,2]),])
cor_err_smooth_total <-cor(pca_refFlat_total$x,pca_smooth_total$x) #columns refer to PCs in pca_smooth$x, while rows refer to PCs in pca_refFlat$x


############################

#plot 3d-plot with MAPD (targets), mrpb (targets) and pc1 for tumors (run code in terminal to get window with movable 3d plot)
save.image("~/.RData")
plot3d(mrpb_targets_refFlat[-as.integer(testset[,2]),2],MAPD_targets_refFlat[-as.integer(testset[,2]),2],pca_refFlat$x[,1],main="Ref flat")
View(cor(pca_refFlat$x[,1:5],as.numeric(MAPD_targets_refFlat[-as.integer(testset[,2]),2])))
View(cor(pca_refFlat$x[,1:5],as.numeric(mrpb_targets_refFlat[-as.integer(testset[,2]),2])))
View(cor(pca_refFlat$x[,1:5],as.numeric(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2])))
View(cor(pca_refFlat$x[,1:5],as.numeric(mrpb_antitargets_refFlat[-as.integer(testset[,2]),2])))


setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/",sep=""))
pdf("3d_MAPD_mrpb_pc1_woY.pdf")
layout(matrix(1,1,1,byrow=TRUE))
scatterplot3d(mrpb_targets_refFlat[-as.integer(testset[,2]),2],MAPD_targets_refFlat[-as.integer(testset[,2]),2],pca_refFlat$x[,1],pch=16,highlight.3d = TRUE)
dev.off()


########################

#Find clusters, classify test sample due to this and calculate error median due to each cluster

#Perform cluster analysis with OPTICS algorithm on PC1-4
x_inter <- pca_refFlat$x[,1:4]
Eps=150
Eps_cl=15
opt <- optics(x_inter, eps=Eps, minPts=5, eps_cl=Eps_cl)
plot(opt)
max(opt$cluster)
#the different clusters
clInd0 <- which(opt$cluster==0)
clInd1 <- which(opt$cluster==1)
clInd2 <- which(opt$cluster==2)
clInd3 <- which(opt$cluster==3)
clInd4 <- which(opt$cluster==4)
clInd5 <- which(opt$cluster==5)
clInd6 <- which(opt$cluster==6)
clInd7 <- which(opt$cluster==7)
clInd8 <- which(opt$cluster==8)
clInd9 <- which(opt$cluster==9)
clInd10 <- which(opt$cluster==10)
layout(1,1,1)
par(mar=c(5,4,4,2)+0.1)
plot(x_inter[,3],x_inter[,2],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points(x_inter[clInd1,3],x_inter[clInd1,2],col="red")
points(x_inter[clInd2,3],x_inter[clInd2,2],col="green")
points(x_inter[clInd3,3],x_inter[clInd3,2],col="blue")
points(x_inter[clInd4,3],x_inter[clInd4,2],col="turquoise")
points(x_inter[clInd5,3],x_inter[clInd5,2],col="purple")
points(x_inter[clInd6,3],x_inter[clInd6,2],col="yellow")
points(x_inter[clInd7,3],x_inter[clInd7,2],col="grey")
points(x_inter[clInd8,3],x_inter[clInd8,2],col="deeppink3")
points(x_inter[clInd9,3],x_inter[clInd9,2],col="darksalmon")
points(x_inter[clInd10,3],x_inter[clInd10,2],col="olivedrab4")
points(x_inter[clInd0,3],x_inter[clInd0,2],col="black")

#plot twice in 3D to see that the clusters look OK
save.image("~/.RData")
plot3d(x_inter[,1],x_inter[,2],x_inter[,3],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(x_inter[clInd1,1],x_inter[clInd1,2],x_inter[clInd1,3],col="red")
points3d(x_inter[clInd2,1],x_inter[clInd2,2],x_inter[clInd2,3],col="green")
points3d(x_inter[clInd3,1],x_inter[clInd3,2],x_inter[clInd3,3],col="blue")
points3d(x_inter[clInd4,1],x_inter[clInd4,2],x_inter[clInd4,3],col="turquoise")
points3d(x_inter[clInd5,1],x_inter[clInd5,2],x_inter[clInd5,3],col="purple")
points3d(x_inter[clInd6,1],x_inter[clInd6,2],x_inter[clInd6,3],col="yellow")
points3d(x_inter[clInd7,1],x_inter[clInd7,2],x_inter[clInd7,3],col="grey")
points3d(x_inter[clInd8,1],x_inter[clInd8,2],x_inter[clInd8,3],col="deeppink3")
points3d(x_inter[clInd9,1],x_inter[clInd9,2],x_inter[clInd9,3],col="darksalmon")
points3d(x_inter[clInd10,1],x_inter[clInd10,2],x_inter[clInd10,3],col="olivedrab4")
points3d(x_inter[clInd0,1],x_inter[clInd0,2],x_inter[clInd0,3],col="black")

plot3d(x_inter[,2],x_inter[,3],x_inter[,4],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(x_inter[clInd1,2],x_inter[clInd1,3],x_inter[clInd1,4],col="red")
points3d(x_inter[clInd2,2],x_inter[clInd2,3],x_inter[clInd2,4],col="green")
points3d(x_inter[clInd3,2],x_inter[clInd3,3],x_inter[clInd3,4],col="blue")
points3d(x_inter[clInd4,2],x_inter[clInd4,3],x_inter[clInd4,4],col="turquoise")
points3d(x_inter[clInd5,2],x_inter[clInd5,3],x_inter[clInd5,4],col="purple")
points3d(x_inter[clInd6,2],x_inter[clInd6,3],x_inter[clInd6,4],col="yellow")
points3d(x_inter[clInd7,2],x_inter[clInd7,3],x_inter[clInd7,4],col="grey")
points3d(x_inter[clInd8,2],x_inter[clInd8,3],x_inter[clInd8,4],col="deeppink3")
points3d(x_inter[clInd9,2],x_inter[clInd9,3],x_inter[clInd9,4],col="darksalmon")
points3d(x_inter[clInd10,2],x_inter[clInd10,3],x_inter[clInd10,4],col="olivedrab4")
points3d(x_inter[clInd0,2],x_inter[clInd0,3],x_inter[clInd0,4],col="black")

#plot each PC in 3d with mrpb and MAPD to see how clusters spread
plot3d(mrpb_targets_refFlat[-as.integer(testset[,2]),2],MAPD_targets_refFlat[-as.integer(testset[,2]),2],x_inter[,1], type="n",main=paste("eps =",Eps," eps_cl =",Eps_cl,"minPts =", 5))
points3d(mrpb_targets_refFlat[-as.integer(testset[,2]),2][clInd1],MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd1],x_inter[clInd1,1],col="red")
points3d(mrpb_targets_refFlat[-as.integer(testset[,2]),2][clInd2],MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd2],x_inter[clInd2,1],col="green")
points3d(mrpb_targets_refFlat[-as.integer(testset[,2]),2][clInd3],MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd3],x_inter[clInd3,1],col="blue")
points3d(mrpb_targets_refFlat[-as.integer(testset[,2]),2][clInd0],MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd0],x_inter[clInd0,1],col="black")

plot3d(mrpb_antitargets_refFlat[-as.integer(testset[,2]),2],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2],x_inter[,1], type="n",main=paste("eps =",Eps," eps_cl =",Eps_cl,"minPts =", 5))
points3d(mrpb_antitargets_refFlat[-as.integer(testset[,2]),2][clInd1],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd1],x_inter[clInd1,1],col="red")
points3d(mrpb_antitargets_refFlat[-as.integer(testset[,2]),2][clInd2],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd2],x_inter[clInd2,1],col="green")
points3d(mrpb_antitargets_refFlat[-as.integer(testset[,2]),2][clInd3],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd3],x_inter[clInd3,1],col="blue")
points3d(mrpb_antitargets_refFlat[-as.integer(testset[,2]),2][clInd0],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd0],x_inter[clInd0,1],col="black")

#check how total read counts and on target rate correlates with clusters
setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/",sep=""))
jpeg("Cluster_sep_totRC_onTargRate.jpeg")
plot(tot_reads[-as.integer(testset[,2]),1],tot_reads[-as.integer(testset[,2]),2],type="n", xlab="total read count",ylab="on target rate",main="Bv6 training set clusters")
points(tot_reads[-as.integer(testset[,2]),1][clInd1],tot_reads[-as.integer(testset[,2]),2][clInd1],col="red")
points(tot_reads[-as.integer(testset[,2]),1][clInd2],tot_reads[-as.integer(testset[,2]),2][clInd2],col="green")
points(tot_reads[-as.integer(testset[,2]),1][clInd3],tot_reads[-as.integer(testset[,2]),2][clInd3],col="blue")
points(tot_reads[-as.integer(testset[,2]),1][clInd0],tot_reads[-as.integer(testset[,2]),2][clInd0],col="black")
legend("topright",c("Cluster 1", "Cluster 2", "Cluster 3", "Outliers"),pch=1, col=c("red","green","blue","black"))
dev.off()

plot(tot_reads[-as.integer(testset[,2]),1],x_inter[,1],type="n")
points(tot_reads[-as.integer(testset[,2]),1][clInd1],x_inter[,1][clInd1],col="red")
points(tot_reads[-as.integer(testset[,2]),1][clInd2],x_inter[,1][clInd2],col="green")
points(tot_reads[-as.integer(testset[,2]),1][clInd3],x_inter[,1][clInd3],col="blue")
points(tot_reads[-as.integer(testset[,2]),1][clInd0],x_inter[,1][clInd0],col="black")

save.image("~/.RData")
plot3d(tot_reads[-as.integer(testset[,2]),1],tot_reads[-as.integer(testset[,2]),2],x_inter[,1],type="n",xlab = "total rc", ylab = "on target rate")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd1],tot_reads[-as.integer(testset[,2]),2][clInd1],x_inter[clInd1,1],col="red")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd2],tot_reads[-as.integer(testset[,2]),2][clInd2],x_inter[clInd2,1],col="green")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd3],tot_reads[-as.integer(testset[,2]),2][clInd3],x_inter[clInd3,1],col="blue")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd0],tot_reads[-as.integer(testset[,2]),2][clInd0],x_inter[clInd0,1],col="black")

plot3d(tot_reads[-as.integer(testset[,2]),1],tot_reads[-as.integer(testset[,2]),2],x_inter[,2],type="n",xlab = "total rc", ylab = "on target rate")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd1],tot_reads[-as.integer(testset[,2]),2][clInd1],x_inter[clInd1,2],col="red")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd2],tot_reads[-as.integer(testset[,2]),2][clInd2],x_inter[clInd2,2],col="green")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd3],tot_reads[-as.integer(testset[,2]),2][clInd3],x_inter[clInd3,2],col="blue")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd0],tot_reads[-as.integer(testset[,2]),2][clInd0],x_inter[clInd0,2],col="black")

plot3d(tot_reads[-as.integer(testset[,2]),1],tot_reads[-as.integer(testset[,2]),2],x_inter[,3],type="n",xlab = "total rc", ylab = "on target rate")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd1],tot_reads[-as.integer(testset[,2]),2][clInd1],x_inter[clInd1,3],col="red")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd2],tot_reads[-as.integer(testset[,2]),2][clInd2],x_inter[clInd2,3],col="green")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd3],tot_reads[-as.integer(testset[,2]),2][clInd3],x_inter[clInd3,3],col="blue")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd0],tot_reads[-as.integer(testset[,2]),2][clInd0],x_inter[clInd0,3],col="black")

plot3d(tot_reads[-as.integer(testset[,2]),1],tot_reads[-as.integer(testset[,2]),2],x_inter[,4],type="n",xlab = "total rc", ylab = "on target rate")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd1],tot_reads[-as.integer(testset[,2]),2][clInd1],x_inter[clInd1,4],col="red")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd2],tot_reads[-as.integer(testset[,2]),2][clInd2],x_inter[clInd2,4],col="green")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd3],tot_reads[-as.integer(testset[,2]),2][clInd3],x_inter[clInd3,4],col="blue")
points3d(tot_reads[-as.integer(testset[,2]),1][clInd0],tot_reads[-as.integer(testset[,2]),2][clInd0],x_inter[clInd0,4],col="black")

View(cor(x_inter,tot_reads[-as.integer(testset[,2]),1]))
View(cor(x_inter,tot_reads[-as.integer(testset[,2]),2]))

ttest_rc_12 <- t.test(tot_reads[-as.integer(testset[,2]),1][clInd1],tot_reads[-as.integer(testset[,2]),1][clInd2])
ttest_rc_13 <- t.test(tot_reads[-as.integer(testset[,2]),1][clInd1],tot_reads[-as.integer(testset[,2]),1][clInd3])

ttest_otr_12 <- t.test(tot_reads[-as.integer(testset[,2]),2][clInd1],tot_reads[-as.integer(testset[,2]),2][clInd2],alternative = "g")
ttest_otr_13 <- t.test(tot_reads[-as.integer(testset[,2]),2][clInd1],tot_reads[-as.integer(testset[,2]),2][clInd3],alternative = "g")

#########

#take a test sample, take error for each bin from err_mat, calc PC comps and classify it due to the given clusters
#first I use a sample which was also used for building the model 
#error for each bin in the test sample
nonadj_err <- err_mat_refFlat[-as.integer(testset[,2]),][1,]
nonadj_err <- nonadj_err[-Ybins]

#the test sample in PC coordinates
pc_test <- (nonadj_err-pca_refFlat$center) %*% pca_refFlat$rotation 
identical(as.vector(pc_test),as.vector(pca_refFlat$x[1,])) #they should be identical since the "test sample" here is the same ad the first in err_mat_refFlat
#classify test sample into a suitable cluster (k nearest neighbor, not weighted)
assign_clust <- knn(train=pca_refFlat$x[,1:4], test=pc_test[1:4], cl=opt$cluster, k=3)
identical(as.integer(levels(assign_clust)[assign_clust]),opt$cluster[1]) #should be identical

#calculate median error of each bin and cluster, adjust bin-log2 of test sample 
cl_samp <- which(opt$cluster==as.integer(levels(assign_clust)[assign_clust]))
med_err <- apply(err_mat_refFlat[cl_samp,-Ybins],2,median)
var_err <- apply(err_mat_refFlat[cl_samp,-Ybins],2,var)
adj_log2 <- nonadj_err - med_err + smooth_mat_refFlat[-as.integer(testset[,2]),-Ybins][1,] #adj_log2=log2-med_err=log2-smooth+smooth-med_err=err-med_err+smooth
adj_err <- adj_log2 - smooth_mat_refFlat[-as.integer(testset[,2]),-Ybins][1,] 

#plot all bins
layout(matrix(1,1,1))
plot(nonadj_err,adj_err,xlim=c(-7,7),ylim=c(-7,7),xlab=paste("non-adjusted err, var=",signif(var(nonadj_err),2),sep=""),ylab=paste("adjusted err, var=", signif(var(adj_err),2),sep=""),main="all bins")
abline(0,1,lty=2)
abline(v=0,h=0,lty=3)
err_diff <- abs(nonadj_err) - abs(adj_err)
length(which(err_diff>0))
length(which(err_diff>0))/length(err_diff) #part of bins where the absolute error is decreased through the transformation
#plot targets
nonadj_err_targ <- nonadj_err[which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="T")]
adj_err_targ <- adj_err[which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="T")]
layout(matrix(1,1,1))
plot(nonadj_err_targ,adj_err_targ,xlim=c(-7,7),ylim=c(-7,7),xlab=paste("non-adjusted err, var=",signif(var(nonadj_err_targ),2),sep=""),ylab=paste("adjusted err, var=", signif(var(adj_err_targ),2),sep=""),main="target bins")
abline(0,1,lty=2)
abline(v=0,h=0,lty=3)
err_diff_targ <- abs(nonadj_err_targ)-abs(adj_err_targ)
length(which(err_diff_targ>0))
length(which(err_diff_targ>0))/length(err_diff_targ)
#plot antitargets
nonadj_err_antitarg <- nonadj_err[which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="AT")]
adj_err_antitarg <- adj_err[which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="AT")]
layout(matrix(1,1,1))
plot(nonadj_err_antitarg,adj_err_antitarg,xlim=c(-7,7),ylim=c(-7,7),xlab=paste("non-adjusted err, var=",signif(var(nonadj_err_antitarg),2),sep=""),ylab=paste("adjusted err, var=", signif(var(adj_err_antitarg),2),sep=""),main="antitarget bins")
abline(0,1,lty=2)
abline(v=0,h=0,lty=3)
err_diff_antitarg <- abs(nonadj_err_antitarg)-abs(adj_err_antitarg)
length(which(err_diff_antitarg>0))
length(which(err_diff_antitarg>0))/length(err_diff_antitarg)

####

#take several test samples (i.e. which have not been involved in building the model)
samples_test <- testset[,1]
nonadj_err_all <- err_mat_refFlat[as.integer(testset[,2]),-Ybins]
centered_nonadj_err_all <- nonadj_err_all-matrix(pca_refFlat$center,nrow=length(nonadj_err_all[,1]),ncol=length(nonadj_err_all[1,]),byrow=TRUE)
pc_test_all <- centered_nonadj_err_all %*% pca_refFlat$rotation
nonadj_log2_all <- nonadj_err_all + smooth_mat_refFlat[as.integer(testset[,2]),-Ybins]

#calculate median error in each bin over each cluster
med_err <- matrix(nrow=(max(opt$cluster)+1),ncol=length(pca_refFlat$rotation[,1]))
rownames(med_err) <- paste("cluster",c(1:3,0))
for (i in 1:max(opt$cluster)) {
  cl_samp <- which(opt$cluster==i)
  med_err[i,] <- apply(err_mat_refFlat[-as.integer(testset[,2]),-Ybins][cl_samp,],2,median)
}
cl_samp <- which(opt$cluster==0)
med_err[max(opt$cluster)+1,] <- apply(err_mat_refFlat[-as.integer(testset[,2]),-Ybins][cl_samp,],2,median)

#pre-allocate 
assign_clust_all <- rep(NA,length(samples_test))
adj_err_all <- matrix(NA,length(samples_test),length(pca_refFlat$rotation[,1]))
adj_log2_all <- matrix(NA,length(samples_test),length(pca_refFlat$rotation[,1]))

#calculate adjusted bin values and errors for each test sample
for (i in 1:length(samples_test)) {
  #the current test sample in PC coordinates
  pc_test <- pc_test_all[i,]

  #classify test sample into a suitable cluster (k nearest neighbor, not weighted)
  assign_clust <- knn(train=pca_refFlat$x[,1:4], test=pc_test[1:4], cl=opt$cluster, k=3)
  assign_clust_all[i] <- as.integer(levels(assign_clust)[assign_clust])
  if (as.integer(levels(assign_clust)[assign_clust])==0) {
    med_err_row <- max(opt$cluster) + 1
  } else {
    med_err_row <- as.integer(levels(assign_clust)[assign_clust])
  }
  
  #adjust bin-log2 of test sample according to median errors in that cluster 
  adj_log2_all[i,] <- nonadj_log2_all[i,] - med_err[med_err_row,] 
  adj_err_all[i,] <- nonadj_err_all[i,] - med_err[med_err_row,] #adj_err=adj_bin-smooth=nonadj_bin-med_err-smooth=nonadj_err-med_err
 
  #plot adj_err vs nonadj_err for all bins
  setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/errors/",sep=""))
  filename <- paste(samples_test[i],"_errors.jpeg",sep="")
  jpeg(filename,width=900,height=900)
  layout(matrix(1,1,1))
  plot(nonadj_err_all[i,],adj_err_all[i,],xlim=c(-7,7),ylim=c(-7,7),xlab="non-adjusted error", ylab="adjusted error", main=paste("Bin errors before and after adjustment, sample ", samples_test[i],", cluster ",as.integer(levels(assign_clust)[assign_clust]),sep=""),type="n",cex.main=0.8)
  abline(0,1,lty=2)
  abline(v=0,h=0,lty=3)
  points(nonadj_err_all[i,which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="AT")],adj_err_all[i,which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="AT")],col="#00868B80",cex=0.5,pch=16)
  points(nonadj_err_all[i,which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="T")],adj_err_all[i,which(dat_cnr_tumors_refFlat$targetType[-Ybins]=="T")],col="#8B225280",cex=0.5,pch=16)
  legend("topleft",c("target bins","antitarget bins", "adj err = non-adj err"),pch=c(16,16,NA),lty=c(NA,NA,2),col=c("#00868B80","#8B225280","black"))
  dev.off()
  
  #plot bin values genome wide
  title1=paste("log2CNratio bins before pca adjustment, sample", samples_test[i])
  title2=paste("log2CNratio bins after pca adjustment, sample", samples_test[i])
  filename <- paste(samples_test[i],"_nonadj_adj_bins.jpeg",sep="")
  setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/adjBins_genomeWide/",sep=""))
  jpeg(filename,width=2000,height=900)
  layout(matrix(c(1,2),2,1,byrow=TRUE))
  plot(dat_cnr_tumors_refFlat$startcoord[-Ybins],nonadj_log2_all[i,],pch=16, xaxt="n", col='#00000050', cex=0.3,ylim=c(-5,5), main=title1,xlab="Chromosome #", ylab="log2 CN ratio",cex.main=0.8)
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"))
  plot(dat_cnr_tumors_refFlat$startcoord[-Ybins],adj_log2_all[i,],pch=16, xaxt="n",col='#00000050', cex=0.3,ylim=c(-5,5), main=title2,xlab="Chromosome #", ylab="log2 CN ratio",cex.main=0.8)
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"))
  dev.off()
}

#code for 3d plot (run in terminal R)
save.image("~/.RData")
plot3d(x_inter[,1],x_inter[,2],x_inter[,3],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(x_inter[clInd1,1],x_inter[clInd1,2],x_inter[clInd1,3],col="red")
points3d(x_inter[clInd2,1],x_inter[clInd2,2],x_inter[clInd2,3],col="green")
points3d(x_inter[clInd3,1],x_inter[clInd3,2],x_inter[clInd3,3],col="blue")
points3d(x_inter[clInd4,1],x_inter[clInd4,2],x_inter[clInd4,3],col="turquoise")
points3d(x_inter[clInd5,1],x_inter[clInd5,2],x_inter[clInd5,3],col="purple")
points3d(x_inter[clInd6,1],x_inter[clInd6,2],x_inter[clInd6,3],col="yellow")
points3d(x_inter[clInd7,1],x_inter[clInd7,2],x_inter[clInd7,3],col="grey")
points3d(x_inter[clInd0,1],x_inter[clInd0,2],x_inter[clInd0,3],col="black")
points3d(pc_test_all[which(assign_clust_all==1),1],pc_test_all[which(assign_clust_all==1),2],pc_test_all[which(assign_clust_all==1),3], size=6,col="red")
points3d(pc_test_all[which(assign_clust_all==2),1],pc_test_all[which(assign_clust_all==2),2],pc_test_all[which(assign_clust_all==2),3], size=6,col="green")
points3d(pc_test_all[which(assign_clust_all==3),1],pc_test_all[which(assign_clust_all==3),2],pc_test_all[which(assign_clust_all==3),3], size=6,col="blue")
points3d(pc_test_all[which(assign_clust_all==0),1],pc_test_all[which(assign_clust_all==0),2],pc_test_all[which(assign_clust_all==0),3], size=6,col="black")

plot3d(x_inter[,3],x_inter[,4],x_inter[,2],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(x_inter[clInd1,3],x_inter[clInd1,4],x_inter[clInd1,2],col="red")
points3d(x_inter[clInd2,3],x_inter[clInd2,4],x_inter[clInd2,2],col="green")
points3d(x_inter[clInd3,3],x_inter[clInd3,4],x_inter[clInd3,2],col="blue")
points3d(x_inter[clInd4,3],x_inter[clInd4,4],x_inter[clInd4,2],col="turquoise")
points3d(x_inter[clInd5,3],x_inter[clInd5,4],x_inter[clInd5,2],col="purple")
points3d(x_inter[clInd6,3],x_inter[clInd6,4],x_inter[clInd6,2],col="yellow")
points3d(x_inter[clInd7,3],x_inter[clInd7,4],x_inter[clInd7,2],col="grey")
points3d(x_inter[clInd0,3],x_inter[clInd0,4],x_inter[clInd0,2],col="black")
points3d(pc_test_all[which(assign_clust_all==1),3],pc_test_all[which(assign_clust_all==1),4],pc_test_all[which(assign_clust_all==1),2], size=6,col="red")
points3d(pc_test_all[which(assign_clust_all==2),3],pc_test_all[which(assign_clust_all==2),4],pc_test_all[which(assign_clust_all==2),2], size=6,col="green")
points3d(pc_test_all[which(assign_clust_all==3),3],pc_test_all[which(assign_clust_all==3),4],pc_test_all[which(assign_clust_all==3),2], size=6,col="blue")
points3d(pc_test_all[which(assign_clust_all==0),3],pc_test_all[which(assign_clust_all==0),4],pc_test_all[which(assign_clust_all==0),2], size=6,col="black")


#save plots data from cluster assignement and error adjustment
setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/",sep=""))
jpeg("cluster_hist.jpeg",width=1000,height=500)
layout(matrix(c(1,2),1,2))
hist(opt$cluster,freq=FALSE,main="Histogram over clusters for training set samples",xlab="cluster",cex.main=0.9)
hist(assign_clust_all,freq=FALSE,main="Histogram over clusters assigned for test set samples",xlab="assigned cluster",cex.main=0.9)
dev.off()

write.table(adj_log2_all,paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/adj_log2_testset.txt",sep=""),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(assign_clust_all,paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/assigend_clust_testset.txt",sep=""),quote=FALSE,row.names=FALSE,col.names = FALSE)
save.image(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/.RData",sep=""))

#####################################################

#Check how the clusters from error-pca spread in the smooth-pca:
plot3d(pca_smooth$x[,1],pca_smooth$x[,2],pca_smooth$x[,3],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(pca_smooth$x[clInd1,1],pca_smooth$x[clInd1,2],pca_smooth$x[clInd1,3],col="red")
points3d(pca_smooth$x[clInd2,1],pca_smooth$x[clInd2,2],pca_smooth$x[clInd2,3],col="green")
points3d(pca_smooth$x[clInd3,1],pca_smooth$x[clInd3,2],pca_smooth$x[clInd3,3],col="blue")
points3d(pca_smooth$x[clInd4,1],pca_smooth$x[clInd4,2],pca_smooth$x[clInd4,3],col="turquoise")
points3d(pca_smooth$x[clInd5,1],pca_smooth$x[clInd5,2],pca_smooth$x[clInd5,3],col="purple")
points3d(pca_smooth$x[clInd6,1],pca_smooth$x[clInd6,2],pca_smooth$x[clInd6,3],col="yellow")
points3d(pca_smooth$x[clInd7,1],pca_smooth$x[clInd7,2],pca_smooth$x[clInd7,3],col="grey")
points3d(pca_smooth$x[clInd8,1],pca_smooth$x[clInd8,2],pca_smooth$x[clInd8,3],col="deeppink3")
points3d(pca_smooth$x[clInd9,1],pca_smooth$x[clInd9,2],pca_smooth$x[clInd9,3],col="darksalmon")
points3d(pca_smooth$x[clInd10,1],pca_smooth$x[clInd10,2],pca_smooth$x[clInd10,3],col="olivedrab4")
points3d(pca_smooth$x[clInd0,1],pca_smooth$x[clInd0,2],pca_smooth$x[clInd0,3],col="black")

plot3d(pca_smooth$x[,6],pca_smooth$x[,4],pca_smooth$x[,5],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(pca_smooth$x[clInd1,6],pca_smooth$x[clInd1,4],pca_smooth$x[clInd1,5],col="red")
points3d(pca_smooth$x[clInd2,6],pca_smooth$x[clInd2,4],pca_smooth$x[clInd2,5],col="green")
points3d(pca_smooth$x[clInd3,6],pca_smooth$x[clInd3,4],pca_smooth$x[clInd3,5],col="blue")
points3d(pca_smooth$x[clInd4,6],pca_smooth$x[clInd4,4],pca_smooth$x[clInd4,5],col="turquoise")
points3d(pca_smooth$x[clInd5,6],pca_smooth$x[clInd5,4],pca_smooth$x[clInd5,5],col="purple")
points3d(pca_smooth$x[clInd6,6],pca_smooth$x[clInd6,4],pca_smooth$x[clInd6,5],col="yellow")
points3d(pca_smooth$x[clInd7,6],pca_smooth$x[clInd7,4],pca_smooth$x[clInd7,5],col="grey")
points3d(pca_smooth$x[clInd8,6],pca_smooth$x[clInd8,4],pca_smooth$x[clInd8,5],col="deeppink3")
points3d(pca_smooth$x[clInd9,6],pca_smooth$x[clInd9,4],pca_smooth$x[clInd9,5],col="darksalmon")
points3d(pca_smooth$x[clInd10,6],pca_smooth$x[clInd10,4],pca_smooth$x[clInd10,5],col="olivedrab4")
points3d(pca_smooth$x[clInd0,6],pca_smooth$x[clInd0,4],pca_smooth$x[clInd0,5],col="black")

###############

#see how mean seg log2CN varies over the genome for each cluster 
mean_0 <- apply(smooth_mat_refFlat[-as.integer(testset[,2]),][clInd0,], 2, mean)
mean_1 <- apply(smooth_mat_refFlat[-as.integer(testset[,2]),][clInd1,], 2, mean)
mean_2 <- apply(smooth_mat_refFlat[-as.integer(testset[,2]),][clInd2,], 2, mean)
mean_3 <- apply(smooth_mat_refFlat[-as.integer(testset[,2]),][clInd3,], 2, mean)
mean_4 <- apply(smooth_mat_refFlat[-as.integer(testset[,2]),][clInd4,], 2, mean)
mean_5 <- apply(smooth_mat_refFlat[-as.integer(testset[,2]),][clInd5,], 2, mean)
mean_6 <- apply(smooth_mat_refFlat[-as.integer(testset[,2]),][clInd6,], 2, mean)
layout(1,1,1)
plot(dat_cnr_tumors_refFlat$startcoord,mean_0,pch=16, xlab="genome pos", ylab="mean seg log2CN",cex=0.6)
points(dat_cnr_tumors_refFlat$startcoord,mean_1,pch=16,cex=0.6, col="red")
points(dat_cnr_tumors_refFlat$startcoord,mean_2,pch=16,cex=0.6, col="green")
points(dat_cnr_tumors_refFlat$startcoord,mean_3,pch=16,cex=0.6, col="blue")
points(dat_cnr_tumors_refFlat$startcoord,mean_4,pch=16,cex=0.6, col="turquoise")
lengths <- c(length(clInd0),length(clInd1),length(clInd2),length(clInd3),length(clInd4))
legend("bottomleft",paste("clust",0:4, "#samples:", lengths), pch=16, col=c("black","red","green","blue","turquoise"))

plot(dat_cnr_tumors_refFlat$startcoord,mean_1,pch=16,cex=0.6, col="red",ylim=c(-2,2))
points(dat_cnr_tumors_refFlat$startcoord,mean_2,pch=16,cex=0.6, col="green")
points(dat_cnr_tumors_refFlat$startcoord,mean_3,pch=16,cex=0.6, col="blue")


#check if the loadings (rotation) matrix bias some parts of the genome
plot(pca_refFlat$rotation[,1])
plot(pca_refFlat$rotation[,2])
plot(pca_refFlat$rotation[,3])
plot(pca_refFlat$rotation[,4])
plot(pca_refFlat$rotation[,5])
plot(pca_refFlat$rotation[,6])

plot(pca_smooth$rotation[,1])
setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/",sep=""))
jpeg(filename="rotation1_genomeWide.jpeg",width=2000,height=900)
medY <- median(pca_refFlat$rotation[Ybins,1])
medNotY <- median(pca_refFlat$rotation[-Ybins,1])
pval <- 
plot(pca_refFlat$rotation[,1],main="Loading matrix over the genome",sub=paste("Median in Y is:",signif(medY,digits=3), "while median in the rest is:", signif(medNotY,digits=3)),ylim=c(-0.3,0.3))
dev.off()
jpeg(filename="rotation2_genomeWide.jpeg",width=2000,height=900)
medY <- median(pca_refFlat$rotation[Ybins,2])
medNotY <- median(pca_refFlat$rotation[-Ybins,2])
plot(pca_refFlat$rotation[,2],main="Loading matrix over the genome",sub=paste("Median in Y is:",signif(medY,digits=3), "while median in the rest is:", signif(medNotY,digits=3)),ylim=c(-0.3,0.3))
dev.off()
jpeg(filename="rotation3_genomeWide.jpeg",width=2000,height=900)
medY <- median(pca_refFlat$rotation[Ybins,3])
medNotY <- median(pca_refFlat$rotation[-Ybins,3])
plot(pca_refFlat$rotation[,3],main="Loading matrix over the genome",sub=paste("Median in Y is:",signif(medY,digits=3), "while median in the rest is:", signif(medNotY,digits=3)),ylim=c(-0.3,0.3))
dev.off()
jpeg(filename="rotation4_genomeWide.jpeg",width=2000,height=900)
medY <- median(pca_refFlat$rotation[Ybins,4])
medNotY <- median(pca_refFlat$rotation[-Ybins,4])
plot(pca_refFlat$rotation[,4],main="Loading matrix over the genome",sub=paste("Median in Y is:",signif(medY,digits=3), "while median in the rest is:", signif(medNotY,digits=3)),ylim=c(-0.3,0.3))
dev.off()
jpeg(filename="rotation5_genomeWide.jpeg",width=2000,height=900)
medY <- median(pca_refFlat$rotation[Ybins,5])
medNotY <- median(pca_refFlat$rotation[-Ybins,5])
plot(pca_refFlat$rotation[,5],main="Loading matrix over the genome",sub=paste("Median in Y is:",signif(medY,digits=3), "while median in the rest is:", signif(medNotY,digits=3)),ylim=c(-0.3,0.3))
dev.off()


######################################

#Instead of classifying samples into clusters and correct for median in each cluster, 
# fit error as 4-dimensional linear function of PC1-4 in each bin and adjust original bin values with this
coeff_fit <- matrix(nrow=5, ncol=length(pca_refFlat$rotation[,1]))
rownames(coeff_fit) <- c("interc","PC1","PC2","PC3","PC4")
#fit function
for (i in 1:length(coeff_fit[1,])) {
  fit_error <- rlm(err_mat_refFlat[-as.integer(testset[,2]),i] ~ x_inter[,1]+x_inter[,2]+x_inter[,3]+x_inter[,4])
  coeff_fit[,i] <- fit_error$coefficients
}
#calc value of error func in each bin for new samples and adjust bin log2 values from this
function_adjustment <- matrix(coeff_fit[1,],nrow=length(testset[,2]), ncol=length(coeff_fit[1,]), byrow=TRUE) + pc_test_all[,1:4] %*% coeff_fit[2:5,]
func_adj_log2_all <- nonadj_log2_all - function_adjustment
#save the adjusted bin values
write.table(func_adj_log2_all,paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/func_adj_log2_testset.txt",sep=""),quote=FALSE,row.names=FALSE,col.names = FALSE)

#save the workspace
save.image(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/.RData",sep=""))
