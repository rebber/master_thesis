#Plot PCA components, MAPD and mrpd for HOPKINS 1000G data, treated in CNVkit with a 50 normals reference
#Created by Rebecka Bergström on 2016-05-30

#load the PCA workspace created with this script (if it's already run once)
#load("/home/rebecka.bergstrom/cnvkit/HOPKINS/hopkins_1000G/pca/.RData")


library(tools)
library(stats)
library(matrixStats)
library(MASS)
library(permute)

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/HOPKINS/hopkins_1000G/"
setwd(paste(path_CNVkit,"pca/",sep=""))

#find the files
files_cnr_tumors <- dir(pattern=".txt", paste(path_CNVkit,"cnrRcComb/addedSmooth/",sep=""), recursive = TRUE)
files_cns_tumors <- dir(pattern=".cns", paste(path_CNVkit,"batchResults",sep=""), recursive = TRUE)
samples <- basename(file_path_sans_ext(files_cnr_tumors))

#load variable testset, defining which samples belong to the test set
load(paste(path_CNVkit,"test30",sep=""))

#preallocate error matrix 
example_bins <- read.table(paste(path_CNVkit,"cnrRcComb/addedSmooth/", files_cnr_tumors[1], sep=""), header = TRUE, as.is=TRUE)
err_mat <- matrix(nrow=length(samples),ncol=length(example_bins$chromosome))
rownames(err_mat) <- samples

#preallocate matrix with smoothed values for each bin
smooth_mat <- matrix(nrow=length(samples),ncol=length(example_bins$chromosome))
rownames(smooth_mat) <- samples

#preallocate for coverage (mrpb) and MAPD 
MAPD_total <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(MAPD_total) <- c("sample","tumors tot")
MAPD_targets <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(MAPD_targets) <- c("sample","tumors targ")
MAPD_antitargets <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),)
colnames(MAPD_antitargets) <- c("sample","tumors atarg")
mrpb_total <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(mrpb_total) <- c("sample","tumors tot")
mrpb_targets <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(mrpb_targets) <- c("sample","tumors targ")
mrpb_antitargets <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(mrpb_antitargets) <- c("sample","tumors atarg")

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
for (i in 1:length(files_cnr_tumors)) { 
  ptm2<-proc.time()
  #the data for the tumors
  dat_cnr_tumors <- read.table(paste(path_CNVkit,"cnrRcComb/addedSmooth/", files_cnr_tumors[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_tumors$startcoord <- rep(NA,length(dat_cnr_tumors$chromosome))
  dat_cnr_tumors$endcoord <- rep(NA,length(dat_cnr_tumors$chromosome))
  
  #add the genome wide coordinates
  for (j in 1:24) {
    index_cnr_tumors <- which(dat_cnr_tumors$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_tumors$startcoord[index_cnr_tumors] <- chr_coord$start[j] + dat_cnr_tumors$start[index_cnr_tumors] - 1 #add the "genome coordinates" for each bin
    dat_cnr_tumors$endcoord[index_cnr_tumors] <- chr_coord$start[j] + dat_cnr_tumors$end[index_cnr_tumors] -1 #add the "genome coordinates" for each bin
  }
  
  #add the total read count and target read part to the tot_reads matrix
  tot_reads[i,1] <- sum(dat_cnr_tumors$reads)
  tot_reads[i,2] <- sum(dat_cnr_tumors$reads[which(dat_cnr_tumors$targetType=="T")])/tot_reads[i,1]
  
  #add the error (log2bin-smooth) from the tumors to the error matrix
  err_mat[i,] <- dat_cnr_tumors$log2 - dat_cnr_tumors$smooth
  
  #add the smooth values for each bin to the smooth matrix
  smooth_mat[i,] <- dat_cnr_tumors$smooth
  
  #find median read depths per bin (mrpb) and store
  mrpb_total[i,2] <- median(dat_cnr_tumors$reads)
  mrpb_targets[i,2] <- median(dat_cnr_tumors$reads[which(dat_cnr_tumors$targetType=="T")])
  mrpb_antitargets[i,2] <- median(dat_cnr_tumors$reads[which(dat_cnr_tumors$targetType=="AT")])
  print(paste("The median read depth per bin in sample", samples[i], "for tumors is", mrpb_targets[i,2], "reads per bin in targets,", mrpb_antitargets[i,2], "reads per bin in antitargets and",mrpb_total[i,2], "reads per bin in total."))
  
  #find the median absolute pairwise difference (MAPD) for the log2CN-ratios and store
  MAPD_total[i,2] <- median(abs(diff(dat_cnr_tumors$log2)))
  MAPD_targets[i,2] <- median(abs(diff(dat_cnr_tumors$log2[which(dat_cnr_tumors$targetType=="T")])))
  MAPD_antitargets[i,2] <- median(abs(diff(dat_cnr_tumors$log2[which(dat_cnr_tumors$targetType=="AT")])))
  
  
  proc.time()-ptm2
}
proc.time()-ptm
#end of for-loop for the 188 samples


#count bins in PTEN:
length(which(dat_cnr_tumors$startcoord>PTEN_genomestart & dat_cnr_tumors$startcoord<PTEN_genomeend & dat_cnr_tumors$targetType=="T"))

#tables of coverage
mrpb_summary <- cbind(mrpb_total,mrpb_targets[,2],mrpb_antitargets[,2],signif(as.numeric(mrpb_targets[,2])/as.numeric(mrpb_antitargets[,2]),digits=3))
colnames(mrpb_summary) <- c("sample","tot", "targets", "antitargets", "ratio_t/at")
mrpb_summary <- rbind(mrpb_summary, c("medians",as.character(apply(matrix(as.double(mrpb_summary[,2:5]),length(mrpb_summary[,1]),4),2,median))))
#write.table(mrpb_summary,file=paste(path_CNVkit,"Rplots/smooth/randTestSet/MAPD/mrpb_pca.txt",sep=""),quote=FALSE,row.names=FALSE)

MAPD_summary <- cbind(MAPD_total,MAPD_targets[,2],MAPD_antitargets[,2],signif(as.numeric(MAPD_targets[,2])/as.numeric(MAPD_antitargets[,2]),digits=3))
colnames(MAPD_summary) <- c("sample","tot", "targets", "antitargets", "ratio_t/at")
MAPD_summary <- rbind(MAPD_summary, c("medians",as.character(apply(matrix(as.double(MAPD_summary[,2:5]),length(MAPD_summary[,1]),4),2,median))))
#write.table(MAPD_summary,file=paste(path_CNVkit,"Rplots/smooth/randTestSet/MAPD/MAPD_pca.txt",sep=""),quote=FALSE,row.names=FALSE)

#######################################

#plot MAPD plots
#setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/MAPD",sep=""))
reads=1:1400
mapd<-c()
for (i in 1:1400) mapd[i]=median(abs(diff(log2(rpois(10000,reads[i])/reads[i]))))
#filename_MAPD <- "MAPD_HOPKINS1000G.pdf"
#pdf(filename_MAPD,width=9,height=9)
layout(matrix(1,1,1,byrow=TRUE))
plot(reads,mapd,xlab='median #reads per bin',ylab='median absolute pairwise difference (log2CNratio)',type='l',main='Minimum stochastic noise, Bv6 all 334 tumors',log="x",ylim=c(0,2))
points(reads, mapd*sqrt(2),type="l",lty=2)
#points(as.double(mrpb_total_ref5[,2]),MAPD_total_ref5[,2],col="red")
#points(as.double(mrpb_total_ref5[,3]),MAPD_total_ref5[,3],col="blue")
cols <- c("red","blue")
points(as.double(mrpb_targets[,2]),MAPD_targets[,2],col=cols[1],pch=16)
points(as.double(mrpb_antitargets[,2]),MAPD_antitargets[,2],col=cols[2],pch=16)
legend("topright", c("targets","antitargets","Poisson distribution","Poisson*sqrt(2)"),col=c(cols,"black","black"),pch=c(rep(16,2),NA,NA),lty=c(rep(NA,2),1,2))
#dev.off()

#########################

#PCA

#remove all bins in Y chromosome since they confuse the PCA, not neccessary for HOPKINS
#Ybins <- which(dat_cnr_tumors$startcoord >= chr_coord$start[24])

#do PCA on the training set (i.e. all samples except the 30 in test set)
#pca <- prcomp(err_mat[-as.integer(testset[,2]),-Ybins])
pca <- prcomp(err_mat[-as.integer(testset[,2]),])

filename<-"PCA_SD_per_PC.jpeg"
jpeg(filename,width=9,height=9,quality=90,units="in",res=900)
layout(matrix(1,1,1))
barplot(pca$sdev)
dev.off()

filename<-"PCvsPC_1-13.jpeg"
jpeg(filename,width=9,height=9,quality=90,units="in",res=900)
layout(matrix(c(1:12),3,4, byrow = TRUE))
par(mar=c(5,4,4,2)+0.1)
plot(pca$x[,1],pca$x[,2],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,2],pca$x[,3],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,3],pca$x[,4],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,4],pca$x[,5],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,5],pca$x[,6],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,6],pca$x[,7],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,7],pca$x[,8],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,8],pca$x[,9],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,9],pca$x[,10],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,10],pca$x[,11],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,11],pca$x[,12],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,12],pca$x[,13],xlim=c(-400,400),ylim=c(-400,400))
dev.off()

filename<-"PCvsPC_13-25.jpeg"
jpeg(filename,width=9,height=9,quality=90,units="in",res=900)
layout(matrix(c(1:12),3,4, byrow = TRUE))
par(mar=c(5,4,4,2)+0.1)
plot(pca$x[,13],pca$x[,14],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,14],pca$x[,15],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,15],pca$x[,16],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,16],pca$x[,17],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,17],pca$x[,18],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,18],pca$x[,19],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,19],pca$x[,20],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,20],pca$x[,21],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,21],pca$x[,22],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,22],pca$x[,23],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,23],pca$x[,24],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,24],pca$x[,25],xlim=c(-400,400),ylim=c(-400,400))
dev.off()

filename<-"PCvsPC_25-37.jpeg"
jpeg(filename,width=9,height=9,quality=90,units="in",res=900)
layout(matrix(c(1:12),3,4, byrow = TRUE))
par(mar=c(5,4,4,2)+0.1)
plot(pca$x[,25],pca$x[,26],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,26],pca$x[,27],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,27],pca$x[,28],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,28],pca$x[,29],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,29],pca$x[,30],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,30],pca$x[,31],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,31],pca$x[,32],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,32],pca$x[,33],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,33],pca$x[,34],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,34],pca$x[,35],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,35],pca$x[,36],xlim=c(-400,400),ylim=c(-400,400))
plot(pca$x[,36],pca$x[,37],xlim=c(-400,400),ylim=c(-400,400))
dev.off()

###########

#plot error for each target bin towards pc 1:8 
#targ <- which(dat_cnr_tumors$targetType[-Ybins]=="T") 
targ <- which(dat_cnr_tumors$targetType=="T") 
err_mat_targ <- err_mat[,targ]
for(i in 1:length(err_mat[1,targ])) {
  layout(matrix(c(1:8),2,4,byrow = TRUE))
  plot(pca$x[,1],err_mat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n",main=paste("target bin ",i))
  text(pca$x[,1],err_mat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca$x[,2],err_mat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca$x[,2],err_mat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca$x[,3],err_mat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca$x[,3],err_mat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca$x[,4],err_mat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca$x[,4],err_mat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca$x[,5],err_mat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca$x[,5],err_mat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca$x[,6],err_mat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca$x[,6],err_mat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca$x[,7],err_mat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca$x[,7],err_mat_targ[-as.integer(testset[,2]),i],label=1:281)
  plot(pca$x[,8],err_mat_targ[-as.integer(testset[,2]),i],ylim=c(-1,1),type="n")
  text(pca$x[,8],err_mat_targ[-as.integer(testset[,2]),i],label=1:281)
  browser()
}


######################

#PCA for "negative control"
#shuffle samples within each bin -> same stdev as before
err_mat_shuff <- err_mat[-as.integer(testset[,2]),]
ptm <- proc.time()
err_mat_shuff <- apply(err_mat_shuff,2,sample)
proc.time()-ptm
#pca_shuff <- prcomp(err_mat_shuff[,-Ybins])
pca_shuff <- prcomp(err_mat_shuff)
layout(matrix(1,1,1))
barplot(pca_shuff$sdev)
layout(matrix(c(1:4),2,2, byrow = TRUE))
plot(pca_shuff$x[,1],pca_shuff$x[,2],xlab="PC1",ylab="PC2",main="PCA, HOPKINS1000G samples")
plot(pca_shuff$x[,2],pca_shuff$x[,3],xlab="PC2",ylab="PC3",main="PCA, HOPKINS1000G samples")
plot(pca_shuff$x[,3],pca_shuff$x[,4],xlab="PC3",ylab="PC4",main="PCA, HOPKINS1000G samples")
plot(pca_shuff$x[,4],pca_shuff$x[,5],xlab="PC4",ylab="PC5",main="PCA, HOPKINS1000G samples")
layout(matrix(c(1:4),2,2, byrow = TRUE))
plot(pca_shuff$x[,5],pca_shuff$x[,6],xlab="PC5",ylab="PC6",main="PCA, HOPKINS1000G samples")
plot(pca_shuff$x[,6],pca_shuff$x[,7],xlab="PC6",ylab="PC7",main="PCA, HOPKINS1000G samples")
plot(pca_shuff$x[,7],pca_shuff$x[,8],xlab="PC7",ylab="PC8",main="PCA, HOPKINS1000G samples")
plot(pca_shuff$x[,8],pca_shuff$x[,9],xlab="PC8",ylab="PC9",main="PCA, HOPKINS1000G samples")

barplot(pca$sdev-pca_shuff$sdev)

#compare our data and control
filename<-"compNegCtrl_SDvsSD.jpeg"
jpeg(filename,width=9,height=9,quality=90,units="in",res=900)
layout(1,1,1)
plot(pca_shuff$sdev,pca$sdev,xlab="SD per PC, neg ctrl",ylab="SD per PC, samples",main="SD per PC, samples vs neg ctrl")
fit<-rlm(pca$sdev ~ pca_shuff$sdev)
points(pca_shuff$sdev,fit$fitted.values,type="l")
legend("topleft",c("PCs",paste("Fitted line, slope = ",signif(fit$coefficients[[2]],digits=3),", interc = ",signif(fit$coefficients[[1]],digits=3),sep="")),lty=c(NA,1),pch=c(1,NA),col=1)
dev.off()

filename<-"compNegCtrl_barplots.jpeg"
jpeg(filename,width=9,height=9,quality=90,units="in",res=900)
layout(matrix(c(1,2),1,2,byrow=TRUE))
par(mar=c(5,4,4,2)+0.1)
barplot(pca$sdev,ylim=c(0,110),xlab="PC",ylab="Standard deviation",main="SD per PC, pca on samples")
barplot(pca_shuff$sdev,ylim=c(0,110),xlab="PC",ylab="Standard deviation",main="SD per PC, pca on neg ctrl")
dev.off()

####################

#check if the loadings (rotation) matrix bias some parts of the genome
filename<-"loadings_PC1-12.jpeg"
jpeg(filename,width=9,height=9,quality=90,units="in",res=900)
layout(matrix(1:12,3,4,byrow=TRUE))
plot(pca$rotation[,1],ylim=c(-0.05,0.05))
plot(pca$rotation[,2],ylim=c(-0.05,0.05))
plot(pca$rotation[,3],ylim=c(-0.05,0.05))
plot(pca$rotation[,4],ylim=c(-0.05,0.05))
plot(pca$rotation[,5],ylim=c(-0.05,0.05))
plot(pca$rotation[,6],ylim=c(-0.05,0.05))
plot(pca$rotation[,7],ylim=c(-0.05,0.05))
plot(pca$rotation[,8],ylim=c(-0.05,0.05))
plot(pca$rotation[,9],ylim=c(-0.05,0.05))
plot(pca$rotation[,10],ylim=c(-0.05,0.05))
plot(pca$rotation[,11],ylim=c(-0.05,0.05))
plot(pca$rotation[,12],ylim=c(-0.05,0.05))
dev.off()

filename<-"loadings_PC13-24.jpeg"
jpeg(filename,width=9,height=9,quality=90,units="in",res=900)
layout(matrix(1:12,3,4,byrow=TRUE))
plot(pca$rotation[,13],ylim=c(-0.05,0.05))
plot(pca$rotation[,14],ylim=c(-0.05,0.05))
plot(pca$rotation[,15],ylim=c(-0.05,0.05))
plot(pca$rotation[,16],ylim=c(-0.05,0.05))
plot(pca$rotation[,17],ylim=c(-0.05,0.05))
plot(pca$rotation[,18],ylim=c(-0.05,0.05))
plot(pca$rotation[,19],ylim=c(-0.05,0.05))
plot(pca$rotation[,20],ylim=c(-0.05,0.05))
plot(pca$rotation[,21],ylim=c(-0.05,0.05))
plot(pca$rotation[,22],ylim=c(-0.05,0.05))
plot(pca$rotation[,23],ylim=c(-0.05,0.05))
plot(pca$rotation[,24],ylim=c(-0.05,0.05))
dev.off()

############################

#plot 3d-plot with MAPD (targets), mrpb (targets) and pc1 for tumors (run code in terminal to get window with movable 3d plot)
save.image("~/.RData")
plot3d(mrpb_targets[-as.integer(testset[,2]),2],MAPD_targets[-as.integer(testset[,2]),2],pca$x[,1],main="Ref flat")
View(cor(pca$x[,1:5],as.numeric(MAPD_targets[-as.integer(testset[,2]),2])))
View(cor(pca$x[,1:5],as.numeric(mrpb_targets[-as.integer(testset[,2]),2])))
View(cor(pca$x[,1:5],as.numeric(MAPD_antitargets[-as.integer(testset[,2]),2])))
View(cor(pca$x[,1:5],as.numeric(mrpb_antitargets[-as.integer(testset[,2]),2])))


######################################

#the interesting PCs
x_inter <- pca$x[,1:24]

#take test samples (i.e. which have not been involved in building the model)
samples_test <- testset[,1]
#nonadj_err_all <- err_mat[as.integer(testset[,2]),-Ybins]
nonadj_err_all <- err_mat[as.integer(testset[,2]),]
centered_nonadj_err_all <- nonadj_err_all-matrix(pca$center,nrow=length(nonadj_err_all[,1]),ncol=length(nonadj_err_all[1,]),byrow=TRUE)
pc_test_all <- centered_nonadj_err_all %*% pca$rotation
#nonadj_log2_all <- nonadj_err_all + smooth_mat[as.integer(testset[,2]),-Ybins]
nonadj_log2_all <- nonadj_err_all + smooth_mat[as.integer(testset[,2]),]


#fit error as 24-dimensional linear function of PC1-24 in each bin and adjust original bin values with this
coeff_fit <- matrix(nrow=25, ncol=length(pca$rotation[,1]))
rownames(coeff_fit) <- c("interc",paste("PC",1:24,sep=""))
#coerce formula for function fitting
funct_formula <- paste("err_mat[-as.integer(testset[,2]),i] ~ ",paste("x_inter[,",1:23,"]","+",sep="",collapse=""),"x_inter[,24]",sep="")
#fit function
for (i in 1:length(coeff_fit[1,])) {
  fit_error <- rlm(err_mat[-as.integer(testset[,2]),i] ~ x_inter[,1]+x_inter[,2]+x_inter[,3]+x_inter[,4]+x_inter[,5]+x_inter[,6]+
                     x_inter[,7]+x_inter[,8]+x_inter[,9]+x_inter[,10]+x_inter[,11]+x_inter[,12]+x_inter[,13]+x_inter[,14]+x_inter[,15]+
                     x_inter[,16]+x_inter[,17]+x_inter[,18]+x_inter[,19]+x_inter[,20]+x_inter[,21]+x_inter[,22]+x_inter[,23]+x_inter[,24])
  coeff_fit[,i] <- fit_error$coefficients
}
#calc value of error func in each bin for new samples and adjust bin log2 values from this
function_adjustment <- matrix(coeff_fit[1,],nrow=length(testset[,2]), ncol=length(coeff_fit[1,]), byrow=TRUE) + pc_test_all[,1:24] %*% coeff_fit[2:25,]
func_adj_log2_all <- nonadj_log2_all - function_adjustment
#save the adjusted bin values
write.table(func_adj_log2_all,paste(path_CNVkit,"pca/func_adj_log2_testset.txt",sep=""),quote=FALSE,row.names=FALSE,col.names = FALSE)

#save the workspace
save.image(paste(path_CNVkit,"pca/.RData",sep=""))
