#Calculate MAPD and mrpb for Bv6 test samples with three diff refs: 10 tumors, 5 normals and 327 normals respectively
#Created by Rebecka Bergström on 2016-04-14

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/"
load(paste(path_CNVkit,"test50/test50",sep="")) #load testset 
samples <- testset[,1]

files_cnr_10tum <- paste("ref10tum/cnrRcComb/addedSmooth/",samples,".txt",sep="")
files_cnr_5norm <- paste("ref5norm/cnrRcComb/addedSmooth/",samples,".txt",sep="")
files_cnr_327norm <- paste("ref327norm/cnrRcComb/addedSmooth/",samples,".txt",sep="")

#preallocate MAPD & mrpb matrices
MAPD_10tum <- matrix(nrow=length(samples),ncol=2)
colnames(MAPD_10tum) <- c("targets","antitargets")
rownames(MAPD_10tum) <- samples
MAPD_5norm <- matrix(nrow=length(samples),ncol=2)
colnames(MAPD_5norm) <- c("targets","antitargets")
rownames(MAPD_5norm) <- samples
MAPD_327norm <- matrix(nrow=length(samples),ncol=2)
colnames(MAPD_327norm) <- c("targets","antitargets")
rownames(MAPD_327norm) <- samples
mrpb_10tum <- matrix(nrow=length(samples),ncol=2)
colnames(mrpb_10tum) <- c("targets","antitargets")
rownames(mrpb_10tum) <- samples
mrpb_5norm <- matrix(nrow=length(samples),ncol=2)
colnames(mrpb_5norm) <- c("targets","antitargets")
rownames(mrpb_5norm) <- samples
mrpb_327norm <- matrix(nrow=length(samples),ncol=2)
colnames(mrpb_327norm) <- c("targets","antitargets")
rownames(mrpb_327norm) <- samples

for (i in 1:length(samples)) {
  dat_cnr_10tum <- read.table(paste(path_CNVkit,files_cnr_10tum[i],sep=""),header=TRUE,as.is=TRUE)
  dat_cnr_5norm <- read.table(paste(path_CNVkit,files_cnr_5norm[i],sep=""),header=TRUE,as.is=TRUE)
  dat_cnr_327norm <- read.table(paste(path_CNVkit,files_cnr_327norm[i],sep=""),header=TRUE,as.is=TRUE)
  
  #find median read depths per bin (mrpb) and store
  mrpb_10tum[i,1] <- median(dat_cnr_10tum$reads[which(dat_cnr_10tum$targetType=="T")])
  mrpb_10tum[i,2] <- median(dat_cnr_10tum$reads[which(dat_cnr_10tum$targetType=="AT")])
  mrpb_5norm[i,1] <- median(dat_cnr_5norm$reads[which(dat_cnr_5norm$targetType=="T")])
  mrpb_5norm[i,2] <- median(dat_cnr_5norm$reads[which(dat_cnr_5norm$targetType=="AT")])
  mrpb_327norm[i,1] <- median(dat_cnr_327norm$reads[which(dat_cnr_327norm$targetType=="T")])
  mrpb_327norm[i,2] <- median(dat_cnr_327norm$reads[which(dat_cnr_327norm$targetType=="AT")])
  
  #find the median absolute pairwise difference (MAPD) for the log2CN-ratios and store
  MAPD_10tum[i,1] <- median(abs(diff(dat_cnr_10tum$log2[which(dat_cnr_10tum$targetType=="T")])))
  MAPD_10tum[i,2] <- median(abs(diff(dat_cnr_10tum$log2[which(dat_cnr_10tum$targetType=="AT")])))
  MAPD_5norm[i,1] <- median(abs(diff(dat_cnr_5norm$log2[which(dat_cnr_5norm$targetType=="T")])))
  MAPD_5norm[i,2] <- median(abs(diff(dat_cnr_5norm$log2[which(dat_cnr_5norm$targetType=="AT")])))
  MAPD_327norm[i,1] <- median(abs(diff(dat_cnr_327norm$log2[which(dat_cnr_327norm$targetType=="T")])))
  MAPD_327norm[i,2] <- median(abs(diff(dat_cnr_327norm$log2[which(dat_cnr_327norm$targetType=="AT")])))
  
}

write.table(mrpb_10tum,file=paste(path_CNVkit,"MAPDcomp/mrpb_ref10tum.txt",sep=""),quote=FALSE,sep="\t")
write.table(mrpb_5norm,file=paste(path_CNVkit,"MAPDcomp/mrpb_ref5norm.txt",sep=""),quote=FALSE,sep="\t")
write.table(mrpb_327norm,file=paste(path_CNVkit,"MAPDcomp/mrpb_ref327norm.txt",sep=""),quote=FALSE,sep="\t")
write.table(MAPD_10tum,file=paste(path_CNVkit,"MAPDcomp/MAPD_ref10tum.txt",sep=""),quote=FALSE,sep="\t")
write.table(MAPD_5norm,file=paste(path_CNVkit,"MAPDcomp/MAPD_ref5norm.txt",sep=""),quote=FALSE,sep="\t")
write.table(MAPD_327norm,file=paste(path_CNVkit,"MAPDcomp/MAPD_ref327norm.txt",sep=""),quote=FALSE,sep="\t")




