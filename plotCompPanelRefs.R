#Plot the segments and bins from CNVkit for panel data with blood reference vs with tumor reference for Bv7 samples with homozygos PTEN deletion
library(tools)
path_CNVkit_blref <- "/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/"
path_CNVkit_turef <- "/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/tumRef/full_tumors/"
setwd(paste(path_CNVkit_turef,"Rplots/refComp/",sep=""))
samples <- c("PTENdel1","PTENdel2","PTENdel3","PTENdel4") #exchange "PTENdel1" etc for actual sample IDs for PTENdeleted samples

files_cnr_blref <- dir(pattern="cnr", paste(path_CNVkit_blref,"batchResults",sep=""), recursive = TRUE)
files_cns_blref <- dir(pattern="cns", paste(path_CNVkit_blref,"batchResults",sep=""), recursive = TRUE)
files_rc_targets_blref <- dir(pattern="readcount.targets.txt", paste(path_CNVkit_blref,"readCounts",sep=""),recursive=TRUE)
files_rc_antitargets_blref <- dir(pattern="readcount.antitargets.txt", paste(path_CNVkit_blref,"readCounts",sep=""),recursive=TRUE)
files_cnr_turef <- dir(pattern="cnr", paste(path_CNVkit_turef,"batchResults",sep=""), recursive = TRUE)
files_cns_turef <- dir(pattern="cns", paste(path_CNVkit_turef,"batchResults",sep=""), recursive = TRUE)
files_rc_targets_turef <- dir(pattern="readcount.targets.txt", paste(path_CNVkit_turef,"readCounts",sep=""),recursive=TRUE)
files_rc_antitargets_turef <- dir(pattern="readcount.antitargets.txt", paste(path_CNVkit_turef,"readCounts",sep=""),recursive=TRUE)

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

#preallocate for MAPD plot
MAPD_total <- matrix(c(samples,rep(NA,2*length(samples))),length(samples),3)
colnames(MAPD) <- c("sample","bloodRef","tumorRef")
MAPD_targets <- matrix(c(samples,rep(NA,2*length(samples))),length(samples),3)
colnames(MAPD_targets) <- c("sample","bloodRef","tumorRef")
MAPD_antitargets <- matrix(c(samples,rep(NA,2*length(samples))),length(samples),3)
colnames(MAPD_antitargets) <- c("sample","bloodRef","tumorRef")
mrpb_total <- matrix(c(samples,rep(NA,2*length(samples))),length(samples),3)
colnames(mrpb) <- c("sample","bloodRef","tumorRef")
mrpb_targets <- matrix(c(samples,rep(NA,2*length(samples))),length(samples),3)
colnames(mrpb_targets) <- c("sample","bloodRef","tumorRef")
mrpb_antitargets <- matrix(c(samples,rep(NA,2*length(samples))),length(samples),3)
colnames(mrpb_antitargets) <- c("sample","bloodRef","tumorRef")

#run for all samples
ptm<-proc.time()
for (i in 1:length(files_cnr_blref)) {
  #the data from the blood reference
  dat_cnr_blref <- read.table(paste(path_CNVkit_blref,"batchResults/", files_cnr_blref[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_blref$startcoord <- rep(NA, length(dat_cnr_blref$chromosome)) 
  dat_cnr_blref$endcoord <- rep(NA, length(dat_cnr_blref$chromosome)) 
  dat_cnr_blref$reads <- rep(NA, length(dat_cnr_blref$chromosome))
  dat_cnr_blref$targetType <- rep(NA, length(dat_cnr_blref$chromosome))
  dat_cns_blref <- read.table(paste(path_CNVkit_blref,"batchResults/", files_cns_blref[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cns_blref$startcoord <- rep(NA, length(dat_cns_blref$chromosome)) 
  dat_cns_blref$endcoord <- rep(NA, length(dat_cns_blref$chromosome)) 
  dat_rc_targets_blref <- read.table(paste(path_CNVkit_blref,"readCounts/", files_rc_targets_blref[i], sep=""), as.is=TRUE)
  colnames(dat_rc_targets_blref) <- c("chromosome","start","end","gene","readCount")
  dat_rc_targets_blref <- unique(dat_rc_targets_blref)
  dat_rc_antitargets_blref <- read.table(paste(path_CNVkit_blref,"readCounts/", files_rc_antitargets_blref[i], sep=""), as.is=TRUE)
  colnames(dat_rc_antitargets_blref) <- c("chromosome","start","end","gene","readCount")
  dat_rc_antitargets_blref <- unique(dat_rc_antitargets_blref)
  dat_rc_blref <- rbind(dat_rc_targets_blref,dat_rc_antitargets_blref)
  #add the genome wide coordinates for each bin
  for (j in 1:24) {
    index_cnr_blref <- which(dat_cnr_blref$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_blref$startcoord[index_cnr_blref] <- chr_coord$start[j] + dat_cnr_blref$start[index_cnr_blref] - 1 #add the "genome coordinates" for each bin
    dat_cnr_blref$endcoord[index_cnr_blref] <- chr_coord$start[j] + dat_cnr_blref$end[index_cnr_blref] -1 #add the "genome coordinates" for each bin
    index_cns_blref <- which(dat_cns_blref$chromosome==chr_coord$chromosome[j])
    dat_cns_blref$startcoord[index_cns_blref] <- chr_coord$start[j] + dat_cns_blref$start[index_cns_blref] - 1 #add the "genome coordinates" for each segment
    dat_cns_blref$endcoord[index_cns_blref] <- chr_coord$start[j] + dat_cns_blref$end[index_cns_blref] -1 #add the "genome coordinates" for each segment
  }
  #add the raw read count
  for (j in 1:length(dat_cnr_blref$chromosome)) {
    readsTarget <- dat_rc_targets_blref$readCount[which(dat_cnr_blref$chromosome[j]==dat_rc_targets_blref$chromosome & dat_cnr_blref$start[j]==dat_rc_targets_blref$start)]
    if (length(readsTarget) > 0) {
      dat_cnr_blref$reads[j] <- readsTarget
      dat_cnr_blref$targetType[j] <- "T"
    } else {
      readsAntitarget <- dat_rc_antitargets_blref$readCount[which(dat_cnr_blref$chromosome[j]==dat_rc_antitargets_blref$chromosome & dat_cnr_blref$start[j]==dat_rc_antitargets_blref$start)]
      if (length(readsAntitarget) > 0) {
        dat_cnr_blref$reads[j] <- readsAntitarget
        dat_cnr_blref$targetType[j] <- "AT"
      } else {
        print(paste("Bin", j, "in sample", samples[i], "had no read count or target type assigned"))
      }
    }
  }
  
  #the data from the tumor reference
  dat_cnr_turef <- read.table(paste(path_CNVkit_turef,"batchResults/", files_cnr_turef[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_turef$startcoord <- rep(NA, length(dat_cnr_turef$chromosome)) 
  dat_cnr_turef$endcoord <- rep(NA, length(dat_cnr_turef$chromosome)) 
  dat_cnr_turef$reads <- rep(NA, length(dat_cnr_turef$chromosome))
  dat_cnr_turef$targetType <- rep(NA, length(dat_cnr_turef$chromosome))
  dat_cns_turef <- read.table(paste(path_CNVkit_turef,"batchResults/", files_cns_turef[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cns_turef$startcoord <- rep(NA, length(dat_cns_turef$chromosome)) 
  dat_cns_turef$endcoord <- rep(NA, length(dat_cns_turef$chromosome)) 
  dat_rc_targets_turef <- read.table(paste(path_CNVkit_turef,"readCounts/", files_rc_targets_turef[i], sep=""), as.is=TRUE)
  colnames(dat_rc_targets_turef) <- c("chromosome","start","end","gene","readCount")
  dat_rc_targets_turef <- unique(dat_rc_targets_turef)
  dat_rc_antitargets_turef <- read.table(paste(path_CNVkit_turef,"readCounts/", files_rc_antitargets_turef[i], sep=""), as.is=TRUE)
  colnames(dat_rc_antitargets_turef) <- c("chromosome","start","end","gene","readCount")
  dat_rc_antitargets_turef <- unique(dat_rc_antitargets_turef)
  dat_rc_turef <- rbind(dat_rc_targets_turef,dat_rc_antitargets_turef)
  #add the genome wide coordinates
  for (j in 1:24) {
    index_cnr_turef <- which(dat_cnr_turef$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_turef$startcoord[index_cnr_turef] <- chr_coord$start[j] + dat_cnr_turef$start[index_cnr_turef] - 1 #add the "genome coordinates" for each bin
    dat_cnr_turef$endcoord[index_cnr_turef] <- chr_coord$start[j] + dat_cnr_turef$end[index_cnr_turef] -1 #add the "genome coordinates" for each bin
    index_cns_turef <- which(dat_cns_turef$chromosome==chr_coord$chromosome[j])
    dat_cns_turef$startcoord[index_cns_turef] <- chr_coord$start[j] + dat_cns_turef$start[index_cns_turef] - 1 #add the "genome coordinates" for each segment
    dat_cns_turef$endcoord[index_cns_turef] <- chr_coord$start[j] + dat_cns_turef$end[index_cns_turef] -1 #add the "genome coordinates" for each segment
  }
  #add the raw read count for each bin
  for (j in 1:length(dat_cnr_turef$chromosome)) {
    readsTarget <- dat_rc_targets_turef$readCount[which(dat_cnr_turef$chromosome[j]==dat_rc_targets_turef$chromosome & dat_cnr_turef$start[j]==dat_rc_targets_turef$start)]
    if (length(readsTarget) > 0) {
      dat_cnr_turef$reads[j] <- readsTarget
      dat_cnr_turef$targetType[j] <- "T"
    } else {
      readsAntitarget <- dat_rc_antitargets_turef$readCount[which(dat_cnr_turef$chromosome[j]==dat_rc_antitargets_turef$chromosome & dat_cnr_turef$start[j]==dat_rc_antitargets_turef$start)]
      if (length(readsAntitarget) > 0) {
        dat_cnr_turef$reads[j] <- readsAntitarget
        dat_cnr_turef$targetType[j] <- "AT"
      } else {
        print(paste("Bin", j, "in sample", samples[i], "had no read count or target type assigned"))
      }
    }
  }
  
  #find median read depths per bin (mrpb) and store
  mrpb_total[i,2] <- median(dat_cnr_blref$reads)
  mrpb_total[i,3] <- median(dat_cnr_turef$reads)
  mrpb_targets[i,2] <- median(dat_cnr_blref$reads[which(dat_cnr_blref$targetType=="T")])
  mrpb_targets[i,3] <- median(dat_cnr_turef$reads[which(dat_cnr_turef$targetType=="T")])
  mrpb_antitargets[i,2] <- median(dat_cnr_blref$reads[which(dat_cnr_blref$targetType=="AT")])
  mrpb_antitargets[i,3] <- median(dat_cnr_turef$reads[which(dat_cnr_turef$targetType=="AT")])
  print(paste("The median read depth per bin in sample", samples[i], " with blood reference is", mrpb_targets[i,2], "reads per bin in targets,", mrpb_antitargets[i,2], "reads per bin in antitargets and", mrpb_total[i,2], "reads per bin in total."))
  print(paste("The median read depth per bin in sample", samples[i], " with tumor reference is", mrpb_targets[i,3], "reads per bin in targets,", mrpb_antitargets[i,3], "reads per bin in antitargets and",mrpb_total[i,3], "reads per bin in total."))
  
  #find the median absolute pairwise difference (MAPD) for the log2CN-ratios and store
  MAPD_total[i,2] <- median(abs(diff(dat_cnr_blref$log2)))
  MAPD_total[i,3] <- median(abs(diff(dat_cnr_turef$log2)))
  MAPD_targets[i,2] <- median(abs(diff(dat_cnr_blref$log2[which(dat_cnr_blref$targetType=="T")])))
  MAPD_targets[i,3] <- median(abs(diff(dat_cnr_turef$log2[which(dat_cnr_turef$targetType=="T")])))
  MAPD_antitargets[i,2] <- median(abs(diff(dat_cnr_blref$log2[which(dat_cnr_blref$targetType=="AT")])))
  MAPD_antitargets[i,3] <- median(abs(diff(dat_cnr_turef$log2[which(dat_cnr_turef$targetType=="AT")])))
  
  #plot genome wide comparison between panel data with blood reference and tumor reference
  title1 <- paste("log2 copy number CNVkit blood ref Bv7 ",basename(file_path_sans_ext(files_cnr_blref[i])),sep="")
  title2 <- paste("log2 copy number CNVkit tumor ref Bv7 ",basename(file_path_sans_ext(files_cnr_turef[i])),sep="")
  filename <- paste("Bv7_PTENdel_compPanelRefs_",basename(file_path_sans_ext(files_cnr_turef[i])),".cnvkit",".pdf",sep="")
  
  pdf(filename,width=20,height=9.8)
  layout(matrix(c(1,2),2,1,byrow=TRUE))
  plot(dat_cnr_blref$startcoord, dat_cnr_blref$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title1,xlab="Chromosome #", ylab="log2 CN")
  segments(dat_cns_blref$startcoord, dat_cns_blref$log2, dat_cns_blref$endcoord, dat_cns_blref$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=chr_coord$start, label=c(1:22,"X","Y"))
  
  plot(dat_cnr_turef$startcoord, dat_cnr_turef$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title2,xlab="Chromosome #", ylab="log2 CN")
  segments(dat_cns_turef$startcoord, dat_cns_turef$log2, dat_cns_turef$endcoord, dat_cns_turef$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=chr_coord$start, label=c(1:22,"X","Y"))
  dev.off()
  
  #plot qq-plots
  titleqq <- paste("QQ plot log2 CN CNVkit blood ref vs tumor ref Bv7 ",basename(file_path_sans_ext(files_cnr_blref[i])),sep="")
  filenameqq <- paste("QQ_Bv7_PTENdel_PanelRefs_",basename(file_path_sans_ext(files_cnr_turef[i])),".cnvkit",".pdf",sep="")
  pdf(filenameqq,width=20, height=9.8)
  layout(matrix(c(1,2),1,2,byrow = TRUE))
  qqplot(dat_cnr_blref$log2,dat_cnr_turef$log2,main=titleqq,xlab="log2CN blood ref",ylab="log2CN tumor ref")
  abline(0,1)
  legend("topleft","y=x",lty=1)
  qqplot(dat_cnr_blref$log2,dat_cnr_turef$log2,ylim=c(-3.5,3.5),xlim=c(-3.5,3.5),main="zoom in",xlab="log2CN blood ref",ylab="log2CN tumor ref")
  abline(0,1)
  legend("topleft","y=x",lty=1)
  dev.off()
  
}

proc.time()-ptm

#plot MAPD plots
setwd(paste(path_CNVkit_turef,"Rplots/",sep=""))
reads=1:2000
mapd<-c()
for (i in 1:2000) mapd[i]=median(abs(diff(log2(rpois(10000,reads[i])/reads[i]))))
filename_MAPD <- "MAPD_PTENdel_BlRefTumRef.pdf"
pdf(filename_MAPD)
layout(matrix(1,1,1,byrow=TRUE))
plot(reads,mapd,xlab='median #reads per bin',ylab='median absolute pairwise difference (log2CNratio)',type='l',main='Minimum stochastic noise',log="x")
points(as.double(mrpb_total[,2]),MAPD_total[,2],col="red")
points(as.double(mrpb_total[,3]),MAPD_total[,3],col="blue")
points(as.double(mrpb_targets[,2]),MAPD_targets[,2],col="brown3")
points(as.double(mrpb_targets[,3]),MAPD_targets[,3],col="darkorchid4")
points(as.double(mrpb_antitargets[,2]),MAPD_antitargets[,2],col="darkorange3")
points(as.double(mrpb_antitargets[,3]),MAPD_antitargets[,3],col="aquamarine4")
legend("topright", c("all bins, with blood ref","all bins, with tumor ref","targets, with blood ref","targets, with tumor ref","antitargets, with blood ref","antitargets, with tumor ref","Poisson distribution"),col=c("red","blue","brown","darkorchid4","darkorange3","aquamarine4","black"),pch=c(rep(1,6),NA),lty=c(rep(NA,6),1))
dev.off()


plot(density(dat_cnr_blref$reads),xlim=c(-100,2000))
plot(density(dat_cnr_blref$end-dat_cnr_blref$start))
quantile(dat_cnr_blref$end-dat_cnr_blref$start,probs=seq(0,1,0.1))
shortBins <- which((dat_cnr_blref$end-dat_cnr_blref$start)<50000)
median((dat_cnr_blref$end[shortBins]-dat_cnr_blref$start[shortBins]))
quantile((dat_cnr_blref$end[shortBins]-dat_cnr_blref$start[shortBins]),probs=seq(0,1,0.1))

#layout(1)
#qqnorm(dat_cnr_blref$log2[which(dat_cnr_blref$chromosome!="Y")])
#qqline(dat_cnr_blref$log2)
#qqnorm(dat_cnr_turef$log2)
#qqline(dat_cnr_turef$log2)
#qqplot(dat_cnr_blref$log2,dat_cnr_turef$log2,ylim=c(-3.5,3.5),xlim=c(-3.5,3.5))
#abline(0,1)

#qqplot(dat_cnr_blref$log2[which(dat_cnr_blref$chromosome!="Y")],dat_cns_blref$log2[which(dat_cnr_blref$chromosome!="Y")])
#abline(0,1)



