#Plot PCA components, MAPD and mrpd for Bv6 panel data (tumors), treated in CNVkit with a flat reference (log2=0 everywhere)
#Created by Rebecka Bergström on 2016-03-24

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
files_cnr_tumors_refFlat <- dir(pattern="T_panel_v1.txt", paste(path_CNVkit,"cnrRcComb",sep=""), recursive = TRUE)
files_cns_tumors_refFlat <- dir(pattern="T_panel_v1.cns", paste(path_CNVkit,"batchResults",sep=""), recursive = TRUE)

samples <- basename(file_path_sans_ext(files_cnr_tumors_refFlat))

#preallocate error matrix ref flat
example_bins_refFlat <- read.table(paste(path_CNVkit,"cnrRcComb/", files_cnr_tumors_refFlat[1], sep=""), header = TRUE, as.is=TRUE)
err_mat_refFlat <- matrix(nrow=length(samples),ncol=length(example_bins_refFlat$chromosome))
rownames(err_mat_refFlat) <- c(files_cnr_tumors_refFlat)

#preallocate matrix with segmented values for each bin
segm_mat_refFlat <- matrix(nrow=length(samples),ncol=length(example_bins_refFlat$chromosome))
rownames(segm_mat_refFlat) <- c(files_cnr_tumors_refFlat)

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
  #with flat reference:
  #the data for the tumors_refFlat
  dat_cnr_tumors_refFlat <- read.table(paste(path_CNVkit,"cnrRcComb/", files_cnr_tumors_refFlat[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_tumors_refFlat$segLog2 <- rep(NA,length(dat_cnr_tumors_refFlat$chromosome))
  dat_cnr_tumors_refFlat$startcoord <- rep(NA,length(dat_cnr_tumors_refFlat$chromosome))
  dat_cnr_tumors_refFlat$endcoord <- rep(NA,length(dat_cnr_tumors_refFlat$chromosome))
  dat_cns_tumors_refFlat <- read.table(paste(path_CNVkit,"batchResults/", files_cns_tumors_refFlat[i], sep=""), header = TRUE, as.is=TRUE)
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
  err_mat_refFlat[i,] <- dat_cnr_tumors_refFlat$segLog2 - dat_cnr_tumors_refFlat$log2
  
  #add the segmentation values for each bin to the segm matrix
  segm_mat_refFlat[i,] <- dat_cnr_tumors_refFlat$segLog2
  
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
  setwd(paste(path_CNVkit,"Rplots/segmentsPlots/",sep=""))
  
  #the indeces and coordinates for zoom in on PTEN 
  startZoomIn <- PTEN_genomestart - 50*15000 #genome coordinate of the first bin included in zoom in
  endZoomIn <- PTEN_genomeend + 50*15000 #genome coordinate of last bin included in zoom in
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
  title4 <- paste("log2 copy number CNVkit, ref flat, Bv6 tumor, ",basename(file_path_sans_ext(files_cnr_tumors_refFlat[i])),sep="")
  title_zoomIn <- "zoom in on PTEN"
  filename2 <- paste("Bv6_genome_PTEN_",basename(file_path_sans_ext(gsub("-TD1-CS1-capped","",files_cnr_tumors_refFlat[i]))),".pdf",sep="")
  xticks <- as.numeric(chr_coord$start)+(as.numeric(chr_coord$end)-as.numeric(chr_coord$start))/2
  
  #plot for tumors with ref flat 
  #pdf(filename2,width=20,height=5)
  #layout(matrix(c(1,1,2), 1, 3, byrow = TRUE)) 
  #plot tumors ref flat 
  #plot(dat_cnr_tumors_refFlat$startcoord, dat_cnr_tumors_refFlat$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title4,xlab="Chromosome #", ylab="log2 CN ratio")
  #segments(dat_cns_tumors_refFlat$startcoord, dat_cns_tumors_refFlat$log2, dat_cns_tumors_refFlat$endcoord, dat_cns_tumors_refFlat$log2, lwd=4, col="#CD8500")
  #abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  #axis(1, at=xticks, label=c(1:22,"X","Y"))
  #abline(v=PTEN_genomestart,lty=2)
  #text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  #zoom in on PTEN 
  #plot(dat_cnr_tumors_refFlat$startcoord[zoomIn_cnr_tumors_refFlat], dat_cnr_tumors_refFlat$log2[zoomIn_cnr_tumors_refFlat], xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN ratio")
  #segments(segmStartXcoordZoomIn_tumors_refFlat,segmYcoordZoomIn_tumors_refFlat,segmEndXcoordZoomIn_tumors_refFlat,segmYcoordZoomIn_tumors_refFlat, lwd=4, col="#CD8500")
  #abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  #axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7))
  #text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1)
  #dev.off()
  
}
proc.time()-ptm
#end of for-loop for the 281 tumors


#count bins in PTEN:
length(which(dat_cnr_tumors_refFlat$startcoord>PTEN_genomestart & dat_cnr_tumors_refFlat$startcoord<PTEN_genomeend & dat_cnr_tumors_refFlat$targetType=="T"))

#tables of coverage
View(cbind(mrpb_total_refFlat,mrpb_targets_refFlat[,2],mrpb_antitargets_refFlat[,2],signif(as.numeric(mrpb_targets_refFlat[,2])/as.numeric(mrpb_antitargets_refFlat[,2]),digits=3)))
#colMedians(cbind(mrpb_total_refFlat[,2],mrpb_targets_refFlat[,2],mrpb_antitargets_refFlat[,2],signif(as.numeric(mrpb_targets_refFlat[,2])/as.numeric(mrpb_antitargets_refFlat[,2]),digits=3)))

#######################################

#plot MAPD plots
setwd(paste(path_CNVkit,"Rplots/",sep=""))
reads=1:4500
mapd<-c()
for (i in 1:4500) mapd[i]=median(abs(diff(log2(rpois(10000,reads[i])/reads[i]))))
filename_MAPD <- "MAPD_Bv6_Tum.pdf"
pdf(filename_MAPD,width=9,height=9)
layout(matrix(1,1,1,byrow=TRUE))
plot(reads,mapd,xlab='median #reads per bin',ylab='median absolute pairwise difference (log2CNratio)',type='l',main='Minimum stochastic noise, Bv6',log="x",ylim=c(0,1))
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
#since some bins are not covered by a segment their error is NA, which s not tolerated by prcomp()
#therefore exchange the NAs to 0 (the errors should be centered around 0 and just a few extra bins with this value should not be a problem)
#(another alternative would be to replace the NAs with the respective median value for that certain bin, but I don't do that now)
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
  print(paste("Sample",samples[which(pca_refFlat$x[,2]==out2_refFlat[i])],"is outlier with value:",signif(pca_refFlat$x[which(pca_refFlat$x[,2]==out2_refFlat[i]),2],digits=3)))
}
quantile(pca_refFlat$x[,2])
layout(matrix(1,1,1))
barplot(pca_refFlat$sdev)

layout(matrix(c(1:4),2,2, byrow = TRUE))
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
setwd(paste(path_CNVkit,"Rplots/pca",sep=""))
filename_pca <- "pca_Bv6_allTumors.pdf"
pdf(filename_pca,width=9.8,height=9.8)
layout(matrix(c(1:4),2,2, byrow = TRUE))
plot(pca_refFlat$x[,1],pca_refFlat$x[,2],xlab="PC1",ylab="PC2",main="PCA, all Bv6 tumors, flat ref")
plot(pca_refFlat$x[,2],pca_refFlat$x[,3],xlab="PC2",ylab="PC3",main="PCA, all Bv6 tumors, flat ref")
plot(pca_refFlat$x[,3],pca_refFlat$x[,4],xlab="PC3",ylab="PC4",main="PCA, all Bv6 tumors, flat ref")
plot(pca_refFlat$x[,4],pca_refFlat$x[,5],xlab="PC4",ylab="PC5",main="PCA, all Bv6 tumors, flat ref")
layout(matrix(c(1:4),2,2, byrow = TRUE))
plot(pca_refFlat$x[,5],pca_refFlat$x[,6],xlab="PC5",ylab="PC6",main="PCA, all Bv6 tumors, flat ref")
plot(pca_refFlat$x[,6],pca_refFlat$x[,7],xlab="PC6",ylab="PC7",main="PCA, all Bv6 tumors, flat ref")
plot(pca_refFlat$x[,7],pca_refFlat$x[,8],xlab="PC7",ylab="PC8",main="PCA, all Bv6 tumors, flat ref")
plot(pca_refFlat$x[,8],pca_refFlat$x[,9],xlab="PC8",ylab="PC9",main="PCA, all Bv6 tumors, flat ref")
dev.off()


###########

#plot error for each target bin towards pc 1:8 when no samples are removed
targ <- which(dat_cnr_tumors_refFlat$targetType=="T")
err_mat_refFlat_targ <- err_mat_refFlat[,targ]
for(i in 1:length(err_mat_refFlat[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca_refFlat$x[,1],err_mat_refFlat_targ[,i],ylim=c(-1,1),type="n",main=paste("target bin ",i))
  text(pca_refFlat$x[,1],err_mat_refFlat_targ[,i],label=1:281)
  plot(pca_refFlat$x[,2],err_mat_refFlat_targ[,i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,2],err_mat_refFlat_targ[,i],label=1:281)
  plot(pca_refFlat$x[,3],err_mat_refFlat_targ[,i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,3],err_mat_refFlat_targ[,i],label=1:281)
  plot(pca_refFlat$x[,4],err_mat_refFlat_targ[,i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,4],err_mat_refFlat_targ[,i],label=1:281)
  plot(pca_refFlat$x[,5],err_mat_refFlat_targ[,i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,5],err_mat_refFlat_targ[,i],label=1:281)
  plot(pca_refFlat$x[,6],err_mat_refFlat_targ[,i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,6],err_mat_refFlat_targ[,i],label=1:281)
  plot(pca_refFlat$x[,7],err_mat_refFlat_targ[,i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,7],err_mat_refFlat_targ[,i],label=1:281)
  plot(pca_refFlat$x[,8],err_mat_refFlat_targ[,i],ylim=c(-1,1),type="n")
  text(pca_refFlat$x[,8],err_mat_refFlat_targ[,i],label=1:281)
  browser()
}


######################

#PCA for "negative control"
#totally random noise (drawn from normal distr, mean=0, stdev=1:
err_mat_random <- matrix(rnorm(length(samples)*length(example_bins_refFlat$chromosome)), nrow=length(samples),ncol=length(example_bins_refFlat$chromosome))
rownames(err_mat_random) <- c(files_cnr_tumors_refFlat)
pca_random <- prcomp(err_mat_random)
layout(matrix(1,1,1))
barplot(pca_random$sdev)


#shuffle samples within each bin -> same stdev as before
err_mat_shuff <- matrix(nrow=length(samples),ncol=length(example_bins_refFlat$chromosome))
ptm <- proc.time()
err_mat_shuff <- apply(err_mat_refFlat,2,sample)
proc.time()-ptm
pca_shuff <- prcomp(err_mat_shuff)
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

targ <- which(dat_cnr_tumors_refFlat$targetType=="T")
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
barplot(pca_refFlat$sdev-pca_shuff$sdev)
plot(pca_shuff$sdev,pca_refFlat$sdev)

plot(pca_random$sdev,pca_refFlat$sdev,xlim=c(10,14))
fit2<-rlm(pca_refFlat$sdev ~ pca_random$sdev)
points(pca_random$sdev,fit2$fitted.values,type="l")

#save plots
filename_negctrl <- "pca_compNegCtrl_Bv6.pdf"
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


###############
#fit function: error dep on PC1-5 for each target bin 
for (j in 1:length(err_mat_refFlat_targ[1,])) {
  #with pca on error-matrix
  fit_bin <- rlm(x=pca_refFlat$x[,1:5],y=err_mat_refFlat_targ[,j])
  #with pca on neg. ctrl.
  fit_bin_negctrl <- rlm(x=pca_shuff$x[,1:5],y=err_mat_shuff_targ[,j])
  
  #"theoretical error" for target bin j
  theor_err_all <- c()
  theor_err_negctrl_all <- c()
  true_err_all <- c()
  #calc theor error for each sample
  for (i in 1:length(pca_refFlat$x[,1])) {
    pc_sample <- pca_refFlat$x[i,]
    theor_err <- sum(pc_sample[1:5]*fit_bin$coefficients)
    pc_sample_negctrl <- as.vector(err_mat_refFlat[i,] %*% pca_shuff$rotation)
    theor_err_negctrl <- sum(pc_sample_negctrl[1:5]*fit_bin_negctrl$coefficients)
    true_err <- err_mat_refFlat_targ[i,j]
    theor_err_all <- append(theor_err_all,theor_err)
    theor_err_negctrl_all <- append(theor_err_negctrl_all,theor_err_negctrl)
    true_err_all <- append(true_err_all,true_err)
  }
  layout(matrix(1:3,1,3))
  plot(true_err_all,theor_err_all,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),xlab=paste("true err, var=",signif(var(true_err_all),2),sep=""),ylab=paste("theor err, var=", signif(var(theor_err_all),2),sep=""),main=paste("target bin",j))
  plot(true_err_all,theor_err_negctrl_all,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),xlab=paste("true err, var=",signif(var(true_err_all),2),sep=""),ylab=paste("theor err neg ctrl, var=", signif(var(theor_err_negctrl_all),2),sep=""))
  plot(theor_err_negctrl_all,theor_err_all,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),xlab=paste("theor err neg ctrl, var=", signif(var(theor_err_negctrl_all),2),sep=""),ylab=paste("theor err, var=", signif(var(theor_err_all),2),sep=""))
  browser()
}


######################

#Compare with segmented values
#PCA on segm mat
#since some bins are not covered by a segment their error is NA, which s not tolerated by prcomp()
#therefore exchange the NAs to 0 (the errors should be centered around 0 and just a few extra bins with this value should not be a problem)
#(another alternative would be to replace the NAs with the respective median value for that certain bin, but I don't do that now)
ind_na_segm <- which(is.na(segm_mat_refFlat))
segm_mat_refFlat[ind_na_segm] <- 0
#err_mat_refFlat[ind_na_refFlat] <- NA

pca_segm <- prcomp(segm_mat_refFlat)
layout(matrix(1:8,2,4,byrow=TRUE))
plot(pca_refFlat$x[,1],pca_refFlat$x[,2],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_refFlat$x[,1],pca_refFlat$x[,2])
plot(pca_refFlat$x[,2],pca_refFlat$x[,3],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_refFlat$x[,2],pca_refFlat$x[,3])
plot(pca_refFlat$x[,3],pca_refFlat$x[,4],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_refFlat$x[,3],pca_refFlat$x[,4])
plot(pca_refFlat$x[,4],pca_refFlat$x[,5],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_refFlat$x[,4],pca_refFlat$x[,5])
plot(pca_segm$x[,1],pca_segm$x[,2],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_segm$x[,1],pca_segm$x[,2])
plot(pca_segm$x[,2],pca_segm$x[,3],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_segm$x[,2],pca_segm$x[,3])
plot(pca_segm$x[,3],pca_segm$x[,4],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_segm$x[,3],pca_segm$x[,4])
plot(pca_segm$x[,4],pca_segm$x[,5],type="n",xlim=c(-100,100),ylim=c(-100,100))
text(pca_segm$x[,4],pca_segm$x[,5])

cor_err_segm <-cor(pca_refFlat$x,pca_segm$x) #columns refer to PCs in pca_segm$x, while rows refer to PCs in pca_refFlat$x
quantile(abs(cor_err_segm),probs=seq(0,1,0.1))
quantile(abs(cor_err_segm[1:5,]),probs=seq(0,1,0.1)) #look at PC 1-5 in pca_refFlat$x (the interesting PCs) 
which(cor_err_segm==max(cor_err_segm[1:5,]))
length(which(abs(cor_err_segm)>0.1))
length(which(abs(cor_err_segm[1:5,])>0.1))
length(which(abs(cor_err_segm[1,])>0.1))
length(which(abs(cor_err_segm[2,])>0.1))
which(abs(cor_err_segm[1,])>0.1)
which(abs(cor_err_segm[2,])>0.1)
which(abs(cor_err_segm[3,])>0.1)
which(abs(cor_err_segm[4,])>0.1)
which(abs(cor_err_segm[5,])>0.1)
plot(which(abs(cor_err_segm[1,])>0.1),rep(1,length(which(abs(cor_err_segm[1,])>0.1))),ylim=c(0,6),xlim=c(0,285),ylab="PC for errors",xlab="PCs for segm with cor>0.1")
points(which(abs(cor_err_segm[2,])>0.1),rep(2,length(which(abs(cor_err_segm[2,])>0.1))))
points(which(abs(cor_err_segm[3,])>0.1),rep(3,length(which(abs(cor_err_segm[3,])>0.1))))
points(which(abs(cor_err_segm[4,])>0.1),rep(4,length(which(abs(cor_err_segm[4,])>0.1))))
points(which(abs(cor_err_segm[5,])>0.1),rep(5,length(which(abs(cor_err_segm[5,])>0.1))))

length(which(abs(cor_err_segm)>0.2))
length(which(abs(cor_err_segm[1:5,])>0.2))
which(abs(cor_err_segm[1,])>0.2)
which(abs(cor_err_segm[2,])>0.2)
which(abs(cor_err_segm[3,])>0.2)
which(abs(cor_err_segm[4,])>0.2)
which(abs(cor_err_segm[5,])>0.2)
plot(which(abs(cor_err_segm[1,])>0.2),rep(1,length(which(abs(cor_err_segm[1,])>0.2))),ylim=c(0,6),xlim=c(0,25),ylab="PC for errors",xlab="PCs for segm with cor>0.2")
points(which(abs(cor_err_segm[2,])>0.2),rep(2,length(which(abs(cor_err_segm[2,])>0.2))))
points(which(abs(cor_err_segm[3,])>0.2),rep(3,length(which(abs(cor_err_segm[3,])>0.2))))
points(which(abs(cor_err_segm[4,])>0.2),rep(4,length(which(abs(cor_err_segm[4,])>0.2))))
points(which(abs(cor_err_segm[5,])>0.2),rep(5,length(which(abs(cor_err_segm[5,])>0.2))))

plot(pca_segm$x[,1],pca_refFlat$x[,1])
plot(pca_segm$x[,1],pca_refFlat$x[,2])
plot(pca_segm$x[,1],pca_refFlat$x[,3])
plot(pca_segm$x[,1],pca_refFlat$x[,4])
plot(pca_segm$x[,1],pca_refFlat$x[,5])

which(abs(cor_err_segm[1:5,])>0.2)
cor_err_segm[1:5,]

#plot3d(rep(1:281,281),rep(1:281,rep(281,281)),as.vector(cor_err_segm),xlab="PC error",ylab="PC segm",zlab="cor")


#find which samples belong to which cluster (PC1>0 and PC1<0 respectively), calc mean for each bin, plot 
segmpc1_high <- which(pca_segm$x[,1]>0)
segmpc1_low <- which(pca_segm$x[,1]<0)
segm_mat_high <- segm_mat_refFlat[segmpc1_high,]
segm_mat_low <- segm_mat_refFlat[segmpc1_low,]
mean_high <- apply(segm_mat_high, 2, mean)
mean_low <- apply(segm_mat_low, 2, mean)
layout(matrix(c(1,2),2,1,byrow=TRUE))
plot(dat_cnr_tumors_refFlat$startcoord,mean_high,pch=16, xlab="genome pos", ylab="mean seg log2CN")
points(dat_cnr_tumors_refFlat$startcoord,mean_low,pch=16, col="red")
legend("bottomleft",c("segm PC1 >0", "segm PC1 <0"), pch=16, col=c("black","red"))
plot(dat_cnr_tumors_refFlat$startcoord,mean_high,pch=16, xlab="genome pos", ylab="mean seg log2CN",ylim=c(-2,2))
points(dat_cnr_tumors_refFlat$startcoord,mean_low,pch=16, col="red")



############################

#plot 3d-plot with MAPD (targets), mrpb (targets) and pc1 for tumors (run code in terminal to get window with movable 3d plot)
plot3d(mrpb_targets_refFlat[,2],MAPD_targets_refFlat[,2],pca_refFlat$x[,1],main="Ref flat")

setwd(paste(path_CNVkit,"Rplots",sep=""))
pdf("3d_MAPD_mrpb_pc1.pdf")
layout(matrix(1,1,1,byrow=TRUE))
scatterplot3d(mrpb_targets_refFlat[,2],MAPD_targets_refFlat[,2],pca_refFlat$x[,1],pch=16,highlight.3d = TRUE)
dev.off()


########################

#Find clusters, classify new sample due to this and fit error function due to each cluster

#Perform cluster analysis with OPTICS algorithm on PC1-5

#start with only PC1-2 to try out (the main clustering also lies within these two dimensions)
x12 <- pca_refFlat$x[,1:2]
ptm<-proc.time()
Eps=150
Eps_cl=9
opt <- optics(x12, eps=Eps, minPts=3, eps_cl=Eps_cl)
proc.time()-ptm
plot(opt)
x12 <- cbind(x12, opt$cluster)
colnames(x12) <- c(colnames(x12)[1:2], "cluster")
max(opt$cluster)
clInd0 <- which(opt$cluster==0)
clInd1 <- which(opt$cluster==1)
clInd2 <- which(opt$cluster==2)
clInd3 <- which(opt$cluster==3)
clInd4 <- which(opt$cluster==4)
clInd5 <- which(opt$cluster==5)
plot(x12[,1],x12[,2],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points(x12[clInd1,1],x12[clInd1,2],col="red")
points(x12[clInd2,1],x12[clInd2,2],col="blue")
points(x12[clInd3,1],x12[clInd3,2],col="green")
points(x12[clInd4,1],x12[clInd4,2],col="orange")
points(x12[clInd5,1],x12[clInd5,2],col="purple")
points(x12[clInd0,1],x12[clInd0,2],col="black")

#for PC1-5
x_inter <- pca_refFlat$x[,1:5]
ptm<-proc.time()
Eps=150
Eps_cl=17
opt <- optics(x_inter, eps=Eps, minPts=3, eps_cl=Eps_cl)
proc.time()-ptm
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
plot(x_inter[,1],x_inter[,2],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points(x_inter[clInd1,1],x_inter[clInd1,2],col="red")
points(x_inter[clInd2,1],x_inter[clInd2,2],col="green")
points(x_inter[clInd3,1],x_inter[clInd3,2],col="blue")
points(x_inter[clInd4,1],x_inter[clInd4,2],col="turquoise")
points(x_inter[clInd5,1],x_inter[clInd5,2],col="purple")
points(x_inter[clInd6,1],x_inter[clInd6,2],col="yellow")
points(x_inter[clInd7,1],x_inter[clInd7,2],col="grey")
points(x_inter[clInd0,1],x_inter[clInd0,2],col="black")

#plot twice in 3D to see that the clusters look OK
plot3d(x_inter[,1],x_inter[,2],x_inter[,3],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(x_inter[clInd1,1],x_inter[clInd1,2],x_inter[clInd1,3],col="red")
points3d(x_inter[clInd2,1],x_inter[clInd2,2],x_inter[clInd2,3],col="green")
points3d(x_inter[clInd3,1],x_inter[clInd3,2],x_inter[clInd3,3],col="blue")
points3d(x_inter[clInd4,1],x_inter[clInd4,2],x_inter[clInd4,3],col="turquoise")
points3d(x_inter[clInd5,1],x_inter[clInd5,2],x_inter[clInd5,3],col="purple")
points3d(x_inter[clInd6,1],x_inter[clInd6,2],x_inter[clInd6,3],col="yellow")
points3d(x_inter[clInd7,1],x_inter[clInd7,2],x_inter[clInd7,3],col="grey")
points3d(x_inter[clInd0,1],x_inter[clInd0,2],x_inter[clInd0,3],col="black")

plot3d(x_inter[,3],x_inter[,4],x_inter[,5],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(x_inter[clInd1,3],x_inter[clInd1,4],x_inter[clInd1,5],col="red")
points3d(x_inter[clInd2,3],x_inter[clInd2,4],x_inter[clInd2,5],col="green")
points3d(x_inter[clInd3,3],x_inter[clInd3,4],x_inter[clInd3,5],col="blue")
points3d(x_inter[clInd4,3],x_inter[clInd4,4],x_inter[clInd4,5],col="turquoise")
points3d(x_inter[clInd5,3],x_inter[clInd5,4],x_inter[clInd5,5],col="purple")
points3d(x_inter[clInd6,3],x_inter[clInd6,4],x_inter[clInd6,5],col="yellow")
points3d(x_inter[clInd7,3],x_inter[clInd7,4],x_inter[clInd7,5],col="grey")
points3d(x_inter[clInd0,3],x_inter[clInd0,4],x_inter[clInd0,5],col="black")

#take a new sample, calc error for each bin, calc PC comps and classify it due to the given clusters
#first I use a sample which was also used for building the model 
dat_cnr_new <- read.table(paste(path_CNVkit,"cnrRcComb/", files_cnr_tumors_refFlat[1], sep=""), header = TRUE, as.is=TRUE)
dat_cnr_new$segLog2 <- rep(NA,length(dat_cnr_new$chromosome))
dat_cns_new <- read.table(paste(path_CNVkit,"batchResults/", files_cns_tumors_refFlat[1], sep=""), header = TRUE, as.is=TRUE)
for (j in 1:length(dat_cns_new$chromosome)) {
  ind <- which(dat_cnr_new$chromosome==dat_cns_new$chromosome[j] & dat_cnr_new$start>=dat_cns_new$start[j] & dat_cnr_new$end<=dat_cns_new$end[j])
  dat_cnr_new$segLog2[ind] <- dat_cns_new$log2[j]
}
#error for each bin in the new sample
err_new <- dat_cnr_new$segLog2 - dat_cnr_new$log2 
ind_na_new <- which(is.na(err_new))
err_new[ind_na_new] <- 0
#the new sample in PC coordinates
pc_new <- (err_new-pca_refFlat$center) %*% pca_refFlat$rotation 
identical(as.vector(pc_new),as.vector(pca_refFlat$x[1,])) #they should be identical since the "new sample" here is the same ad the first in err_mat_refFlat
#classify new sample into a suitable cluster (k nearest neighbor, not weighted)
assign_clust <- knn(train=pca_refFlat$x, test=pc_new, cl=opt$cluster, k=3)
identical(as.integer(levels(assign_clust)[assign_clust]),opt$cluster[1]) #should be identical

#fit function: error dep on PC1-5 of the samples in the assigned cluster for each bin 
cl_samp <- which(opt$cluster==as.integer(levels(assign_clust)[assign_clust]))
theor_err_all <- c()
true_err_all <- c()
for (j in 1:length(err_mat_refFlat[1,])) {
  #with pca on error-matrix
  fit_bin <- rlm(x=pca_refFlat$x[cl_samp,1:5],y=err_mat_refFlat[cl_samp,j]) #use centered error instead for fitting?
  
  #"theoretical error" for target bin j
  theor_err <- sum(pc_new[1:5]*fit_bin$coefficients)
  true_err <- err_new[j]
  theor_err_all <- append(theor_err_all,theor_err)
  true_err_all <- append(true_err_all,true_err)
}
#plot all bins
layout(matrix(1,1,1))
plot(true_err_all,theor_err_all,xlim=c(-7,7),ylim=c(-7,7),xlab=paste("true err, var=",signif(var(true_err_all),2),sep=""),ylab=paste("theor err, var=", signif(var(theor_err_all),2),sep=""),main="all bins")
abline(0,1,lty=2)
abline(v=0,h=0,lty=3)
err_diff_all <- abs(true_err_all) - abs(theor_err_all)
length(which(err_diff_all>0))
length(which(err_diff_all>0))/length(err_diff_all) #part of bins where the absolute error is decreased through the transformation
#plot targets
true_err_targ <- true_err_all[which(dat_cnr_tumors_refFlat$targetType=="T")]
theor_err_targ <- theor_err_all[which(dat_cnr_tumors_refFlat$targetType=="T")]
layout(matrix(1,1,1))
plot(true_err_targ,theor_err_targ,xlim=c(-7,7),ylim=c(-7,7),xlab=paste("true err, var=",signif(var(true_err_targ),2),sep=""),ylab=paste("theor err, var=", signif(var(theor_err_targ),2),sep=""),main="target bins")
abline(0,1,lty=2)
abline(v=0,h=0,lty=3)
err_diff_targ <- abs(true_err_targ)-abs(theor_err_targ)
length(which(err_diff_targ>0))
length(which(err_diff_targ>0))/length(err_diff_targ)
#plot antitargets
true_err_antitarg <- true_err_all[which(dat_cnr_tumors_refFlat$targetType=="AT")]
theor_err_antitarg <- theor_err_all[which(dat_cnr_tumors_refFlat$targetType=="AT")]
layout(matrix(1,1,1))
plot(true_err_antitarg,theor_err_antitarg,xlim=c(-7,7),ylim=c(-7,7),xlab=paste("true err, var=",signif(var(true_err_antitarg),2),sep=""),ylab=paste("theor err, var=", signif(var(theor_err_antitarg),2),sep=""),main="antitarget bins")
abline(0,1,lty=2)
abline(v=0,h=0,lty=3)
err_diff_antitarg <- abs(true_err_antitarg)-abs(theor_err_antitarg)
length(which(err_diff_antitarg>0))
length(which(err_diff_antitarg>0))/length(err_diff_antitarg)


#take several actually new samples (i.e. which have not been involved in building the model)
files_cnr_new <- dir(pattern=".cnr",paste(path_CNVkit,"batchResultsDevs",sep=""),recursive=TRUE)
files_cns_new <- dir(pattern=".cns",paste(path_CNVkit,"batchResultsDevs",sep=""),recursive=TRUE)
samples_new <- basename(file_path_sans_ext(files_cnr_new))
pc_new_all <-matrix(NA,length(samples_new),281)
assign_clust_all <- rep(NA,length(samples_new))
for (i in 1:length(samples_new)) {
  dat_cnr_new <- read.table(paste(path_CNVkit,"batchResultsDevs/", files_cnr_new[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_new$segLog2 <- rep(NA,length(dat_cnr_new$chromosome))
  dat_cns_new <- read.table(paste(path_CNVkit,"batchResultsDevs/", files_cns_new[i], sep=""), header = TRUE, as.is=TRUE)
  for (j in 1:length(dat_cns_new$chromosome)) {
    ind <- which(dat_cnr_new$chromosome==dat_cns_new$chromosome[j] & dat_cnr_new$start>=dat_cns_new$start[j] & dat_cnr_new$end<=dat_cns_new$end[j])
    dat_cnr_new$segLog2[ind] <- dat_cns_new$log2[j]
  }
  #error for each bin in the new sample
  err_new <- dat_cnr_new$segLog2 - dat_cnr_new$log2 
  ind_na_new <- which(is.na(err_new))
  err_new[ind_na_new] <- 0
  #the new sample in PC coordinates
  pc_new <- (err_new-pca_refFlat$center) %*% pca_refFlat$rotation 
  pc_new_all[i,] <- pc_new
  #classify new sample into a suitable cluster (k nearest neighbor, not weighted)
  assign_clust <- knn(train=pca_refFlat$x, test=pc_new, cl=opt$cluster, k=3)
  assign_clust_all[i] <- as.integer(levels(assign_clust)[assign_clust])
  
  #fit function: error dep on PC1-5 of the samples in the assigned cluster for each bin 
  cl_samp <- which(opt$cluster==as.integer(levels(assign_clust)[assign_clust]))
  theor_err_all <- c()
  true_err_all <- c()
  for (j in 1:length(err_mat_refFlat[1,])) {
    #with pca on error-matrix
    fit_bin <- rlm(x=pca_refFlat$x[cl_samp,1:5],y=err_mat_refFlat[cl_samp,j]) #use centered error instead, for fitting?
    
    #"theoretical error" for target bin j
    theor_err <- sum(pc_new[1:5]*fit_bin$coefficients)
    true_err <- err_new[j]
    theor_err_all <- append(theor_err_all,theor_err)
    true_err_all <- append(true_err_all,true_err)
  }
  #plot all bins
  setwd(paste(path_CNVkit,"Rplots/pca/errors/",sep=""))
  filename <- paste(samples_new[i],"_allBins_errors.jpeg",sep="")
  jpeg(filename)
  layout(matrix(1,1,1))
  plot(true_err_all,theor_err_all,xlim=c(-7,7),ylim=c(-7,7),xlab=paste("true err, var=",signif(var(true_err_all),2),sep=""),ylab=paste("theor err, var=", signif(var(theor_err_all),2),sep=""),main=paste("all bins, sample ", samples_new[i],", cluster ",as.integer(levels(assign_clust)[assign_clust]),sep=""))
  abline(0,1,lty=2)
  abline(v=0,h=0,lty=3)
  dev.off()
  err_diff_all <- abs(true_err_all) - abs(theor_err_all)
  length(which(err_diff_all>0))
  length(which(err_diff_all>0))/length(err_diff_all) #part of bins where the absolute error is decreased through the transformation
  #plot targets
  filename <- paste(samples_new[i],"_targets_errors.jpeg",sep="")
  jpeg(filename)
  true_err_targ <- true_err_all[which(dat_cnr_tumors_refFlat$targetType=="T")]
  theor_err_targ <- theor_err_all[which(dat_cnr_tumors_refFlat$targetType=="T")]
  layout(matrix(1,1,1))
  plot(true_err_targ,theor_err_targ,xlim=c(-7,7),ylim=c(-7,7),xlab=paste("true err, var=",signif(var(true_err_targ),2),sep=""),ylab=paste("theor err, var=", signif(var(theor_err_targ),2),sep=""),main=paste("target bins, sample ", samples_new[i],", cluster ",as.integer(levels(assign_clust)[assign_clust]),sep=""))   
  abline(0,1,lty=2)
  abline(v=0,h=0,lty=3)
  dev.off()
  err_diff_targ <- abs(true_err_targ)-abs(theor_err_targ)
  length(which(err_diff_targ>0))
  length(which(err_diff_targ>0))/length(err_diff_targ)
  #plot antitargets
  filename <- paste(samples_new[i],"_antitargets_errors.jpeg",sep="")
  jpeg(filename)
  true_err_antitarg <- true_err_all[which(dat_cnr_tumors_refFlat$targetType=="AT")]
  theor_err_antitarg <- theor_err_all[which(dat_cnr_tumors_refFlat$targetType=="AT")]
  layout(matrix(1,1,1))
  plot(true_err_antitarg,theor_err_antitarg,xlim=c(-7,7),ylim=c(-7,7),xlab=paste("true err, var=",signif(var(true_err_antitarg),2),sep=""),ylab=paste("theor err, var=", signif(var(theor_err_antitarg),2),sep=""),main=paste("antitarget bins, sample ", samples_new[i],", cluster ",as.integer(levels(assign_clust)[assign_clust]),sep=""))   
  abline(0,1,lty=2)
  abline(v=0,h=0,lty=3)
  dev.off()
  err_diff_antitarg <- abs(true_err_antitarg)-abs(theor_err_antitarg)
  length(which(err_diff_antitarg>0))
  length(which(err_diff_antitarg>0))/length(err_diff_antitarg)
}

#code for 3d plot (run in terminal R)
setwd("~")
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
points3d(pc_new_all[,1],pc_new_all[,2],pc_new_all[,3], size=6,col="orange")

plot3d(x_inter[,3],x_inter[,4],x_inter[,5],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(x_inter[clInd1,3],x_inter[clInd1,4],x_inter[clInd1,5],col="red")
points3d(x_inter[clInd2,3],x_inter[clInd2,4],x_inter[clInd2,5],col="green")
points3d(x_inter[clInd3,3],x_inter[clInd3,4],x_inter[clInd3,5],col="blue")
points3d(x_inter[clInd4,3],x_inter[clInd4,4],x_inter[clInd4,5],col="turquoise")
points3d(x_inter[clInd5,3],x_inter[clInd5,4],x_inter[clInd5,5],col="purple")
points3d(x_inter[clInd6,3],x_inter[clInd6,4],x_inter[clInd6,5],col="yellow")
points3d(x_inter[clInd7,3],x_inter[clInd7,4],x_inter[clInd7,5],col="grey")
points3d(x_inter[clInd0,3],x_inter[clInd0,4],x_inter[clInd0,5],col="black")
points3d(pc_new_all[,3],pc_new_all[,4],pc_new_all[,5], size=6,col="orange")


#####################################################

#Check how the clusters from error-pca spread in the segement-pca:
plot3d(pca_segm$x[,1],pca_segm$x[,2],pca_segm$x[,3],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(pca_segm$x[clInd1,1],pca_segm$x[clInd1,2],pca_segm$x[clInd1,3],col="red")
points3d(pca_segm$x[clInd2,1],pca_segm$x[clInd2,2],pca_segm$x[clInd2,3],col="green")
points3d(pca_segm$x[clInd3,1],pca_segm$x[clInd3,2],pca_segm$x[clInd3,3],col="blue")
points3d(pca_segm$x[clInd4,1],pca_segm$x[clInd4,2],pca_segm$x[clInd4,3],col="turquoise")
points3d(pca_segm$x[clInd5,1],pca_segm$x[clInd5,2],pca_segm$x[clInd5,3],col="purple")
points3d(pca_segm$x[clInd6,1],pca_segm$x[clInd6,2],pca_segm$x[clInd6,3],col="yellow")
points3d(pca_segm$x[clInd7,1],pca_segm$x[clInd7,2],pca_segm$x[clInd7,3],col="grey")
points3d(pca_segm$x[clInd0,1],pca_segm$x[clInd0,2],pca_segm$x[clInd0,3],col="black")

plot3d(pca_segm$x[,3],pca_segm$x[,4],pca_segm$x[,5],type="n", main=paste("eps =",Eps," eps_cl =",Eps_cl))
points3d(pca_segm$x[clInd1,3],pca_segm$x[clInd1,4],pca_segm$x[clInd1,5],col="red")
points3d(pca_segm$x[clInd2,3],pca_segm$x[clInd2,4],pca_segm$x[clInd2,5],col="green")
points3d(pca_segm$x[clInd3,3],pca_segm$x[clInd3,4],pca_segm$x[clInd3,5],col="blue")
points3d(pca_segm$x[clInd4,3],pca_segm$x[clInd4,4],pca_segm$x[clInd4,5],col="turquoise")
points3d(pca_segm$x[clInd5,3],pca_segm$x[clInd5,4],pca_segm$x[clInd5,5],col="purple")
points3d(pca_segm$x[clInd6,3],pca_segm$x[clInd6,4],pca_segm$x[clInd6,5],col="yellow")
points3d(pca_segm$x[clInd7,3],pca_segm$x[clInd7,4],pca_segm$x[clInd7,5],col="grey")
points3d(pca_segm$x[clInd0,3],pca_segm$x[clInd0,4],pca_segm$x[clInd0,5],col="black")

###############

#see how mean seg log2CN varies over the genome for each cluster 
mean_0 <- apply(segm_mat_refFlat[clInd0,], 2, mean)
mean_1 <- apply(segm_mat_refFlat[clInd1,], 2, mean)
mean_2 <- apply(segm_mat_refFlat[clInd2,], 2, mean)
mean_3 <- apply(segm_mat_refFlat[clInd3,], 2, mean)
mean_4 <- apply(segm_mat_refFlat[clInd4,], 2, mean)
mean_5 <- apply(segm_mat_refFlat[clInd5,], 2, mean)
mean_6 <- apply(segm_mat_refFlat[clInd6,], 2, mean)
layout(matrix(c(1,2),2,1,byrow=TRUE))
plot(dat_cnr_tumors_refFlat$startcoord,mean_0,pch=16, xlab="genome pos", ylab="mean seg log2CN",cex=0.6)
points(dat_cnr_tumors_refFlat$startcoord,mean_1,pch=16,cex=0.6, col="red")
points(dat_cnr_tumors_refFlat$startcoord,mean_2,pch=16,cex=0.6, col="green")
points(dat_cnr_tumors_refFlat$startcoord,mean_3,pch=16,cex=0.6, col="blue")
points(dat_cnr_tumors_refFlat$startcoord,mean_4,pch=16,cex=0.6, col="turquoise")
points(dat_cnr_tumors_refFlat$startcoord,mean_5,pch=16,cex=0.6, col="purple")
points(dat_cnr_tumors_refFlat$startcoord,mean_6,pch=16,cex=0.6, col="yellow")
legend("bottomleft",paste("clust",0:6, "#samples#, "), pch=16, col=c("black","red","green","blue","turquoise","purple","yellow"))
plot(dat_cnr_tumors_refFlat$startcoord,mean_0,pch=16, xlab="genome pos", ylab="mean seg log2CN",cex=0.6,ylim=c(-2,2))
points(dat_cnr_tumors_refFlat$startcoord,mean_1,pch=16,cex=0.6, col="red")
points(dat_cnr_tumors_refFlat$startcoord,mean_2,pch=16,cex=0.6, col="green")
points(dat_cnr_tumors_refFlat$startcoord,mean_3,pch=16,cex=0.6, col="blue")
points(dat_cnr_tumors_refFlat$startcoord,mean_4,pch=16,cex=0.6, col="turquoise")
points(dat_cnr_tumors_refFlat$startcoord,mean_5,pch=16,cex=0.6, col="purple")
points(dat_cnr_tumors_refFlat$startcoord,mean_6,pch=16,cex=0.6, col="yellow")

plot(dat_cnr_tumors_refFlat$startcoord,mean_1,pch=16,cex=0.6, col="red",ylim=c(-2,2))
points(dat_cnr_tumors_refFlat$startcoord,mean_2,pch=16,cex=0.6, col="green")
points(dat_cnr_tumors_refFlat$startcoord,mean_3,pch=16,cex=0.6, col="blue")


#check if the loadings (rotation) matrix bias some parts of the genome
plot(pca_segm$rotation[,1])
setwd(paste(path_CNVkit,"Rplots/pca",sep=""))
jpeg(filename="rotation1_genomeWide.jpeg",width=4800,height=480)
plot(pca_refFlat$rotation[,1])
dev.off()
plot(pca_refFlat$rotation[,2])
plot(pca_refFlat$rotation[,3])
plot(pca_refFlat$rotation[,4])
plot(pca_refFlat$rotation[,5])
plot(pca_refFlat$rotation[,6])



