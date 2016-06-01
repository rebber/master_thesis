#Compare non-adjusted bins and segments (directly from CNVkit) with pca-adjusted bins and segments.
#Created by Rebecka Bergström on 2016-04-14

library(tools)

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/"
load(paste(path_CNVkit,"test50/test50",sep="")) #the variable testset

samples <- file_path_sans_ext(dir(pattern="cnr",path=paste(path_CNVkit,"segmAdjBins/segmResults/",sep=""),recursive = TRUE))

files_cnr_nonadj <- paste("cnrRcComb/addedSmooth/",samples,".txt",sep="")
files_cns_nonadj <- paste("allCnr/",samples,".cns",sep="")
files_cnr_adj <- paste("segmAdjBins/segmResults/",samples,".cnr",sep="")
files_cns_adj <- paste("segmAdjBins/segmResults/",samples,".cns",sep="")

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

#preallocate MAPD & mrpb matrices
MAPD_nonadj <- matrix(nrow=length(samples),ncol=2)
colnames(MAPD_nonadj) <- c("targets","antitargets")
rownames(MAPD_nonadj) <- samples
MAPD_adj <- matrix(nrow=length(samples),ncol=2)
colnames(MAPD_adj) <- c("targets","antitargets")
rownames(MAPD_adj) <- samples
mrpb_nonadj <- matrix(nrow=length(samples),ncol=2)
colnames(mrpb_nonadj) <- c("targets","antitargets")
rownames(mrpb_nonadj) <- samples
mrpb_adj <- matrix(nrow=length(samples),ncol=2)
colnames(mrpb_adj) <- c("targets","antitargets")
rownames(mrpb_adj) <- samples


#for all test samlples read in files and plot
for (i in 1:length(samples)) {
  dat_cnr_nonadj <- read.table(paste(path_CNVkit,files_cnr_nonadj[i],sep=""),header=TRUE,as.is=TRUE)
  dat_cnr_nonadj$startcoord <- rep(NA,length(dat_cnr_nonadj$chromosome))
  dat_cnr_nonadj$endcoord <- rep(NA,length(dat_cnr_nonadj$chromosome))
  dat_cns_nonadj <- read.table(paste(path_CNVkit,files_cns_nonadj[i],sep=""),header=TRUE,as.is=TRUE)
  dat_cns_nonadj$startcoord <- rep(NA,length(dat_cns_nonadj$chromosome))
  dat_cns_nonadj$endcoord <- rep(NA,length(dat_cns_nonadj$chromosome))
  dat_cnr_adj <- read.table(paste(path_CNVkit,files_cnr_adj[i],sep=""),header=TRUE,as.is=TRUE)
  dat_cnr_adj$reaads <- dat_cnr_nonadj$reads
  dat_cnr_adj$targetType <- dat_cnr_nonadj$targetType
  dat_cnr_adj$startcoord <- rep(NA,length(dat_cnr_adj$chromosome))
  dat_cnr_adj$endcoord <- rep(NA,length(dat_cnr_adj$chromosome))
  dat_cns_adj <- read.table(paste(path_CNVkit,files_cns_adj[i],sep=""),header=TRUE,as.is=TRUE)
  dat_cns_adj$startcoord <- rep(NA,length(dat_cns_adj$chromosome))
  dat_cns_adj$endcoord <- rep(NA,length(dat_cns_adj$chromosome))
  
  #add the genome wide coordinates
  for (j in 1:24) {
    #add the "genome coordinates" for each bin
    index_cnr <- which(dat_cnr_nonadj$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_nonadj$startcoord[index_cnr] <- chr_coord$start[j] + dat_cnr_nonadj$start[index_cnr] - 1 
    dat_cnr_nonadj$endcoord[index_cnr] <- chr_coord$start[j] + dat_cnr_nonadj$end[index_cnr] -1 
    dat_cnr_adj$startcoord[index_cnr] <- chr_coord$start[j] + dat_cnr_adj$start[index_cnr] - 1 
    dat_cnr_adj$endcoord[index_cnr] <- chr_coord$start[j] + dat_cnr_adj$end[index_cnr] -1 
    #add the "genome coordinates" for each segment
    index_cns_nonadj <- which(dat_cns_nonadj$chromosome==chr_coord$chromosome[j])
    dat_cns_nonadj$startcoord[index_cns_nonadj] <- chr_coord$start[j] + dat_cns_nonadj$start[index_cns_nonadj] - 1 
    dat_cns_nonadj$endcoord[index_cns_nonadj] <- chr_coord$start[j] + dat_cns_nonadj$end[index_cns_nonadj] -1 
    index_cns_adj <- which(dat_cns_adj$chromosome==chr_coord$chromosome[j])
    dat_cns_adj$startcoord[index_cns_adj] <- chr_coord$start[j] + dat_cns_adj$start[index_cns_adj] - 1 
    dat_cns_adj$endcoord[index_cns_adj] <- chr_coord$start[j] + dat_cns_adj$end[index_cns_adj] -1 
  }
  
  #find median read depths per bin (mrpb) and store
  mrpb_nonadj[i,1] <- median(dat_cnr_nonadj$reads[which(dat_cnr_nonadj$targetType=="T")])
  mrpb_nonadj[i,2] <- median(dat_cnr_nonadj$reads[which(dat_cnr_nonadj$targetType=="AT")])
  mrpb_adj[i,] <- mrpb_nonadj[i,]
  #find the median absolute pairwise difference (MAPD) for the log2CN-ratios and store
  MAPD_nonadj[i,1] <- median(abs(diff(dat_cnr_nonadj$log2[which(dat_cnr_nonadj$targetType=="T")])))
  MAPD_nonadj[i,2] <- median(abs(diff(dat_cnr_nonadj$log2[which(dat_cnr_nonadj$targetType=="AT")])))
  MAPD_adj[i,1] <- median(abs(diff(dat_cnr_adj$log2[which(dat_cnr_adj$targetType=="T")])))
  MAPD_adj[i,2] <- median(abs(diff(dat_cnr_adj$log2[which(dat_cnr_adj$targetType=="AT")])))
  
  
  #plot genomewide and PTEN
  setwd(paste(path_CNVkit,"segmAdjBins/Rplots/compSegm/",sep=""))
  
  #the indeces and coordinates for zoom in on PTEN 
  startZoomIn <- PTEN_genomestart - 50*15000 #genome coordinate of the first bin included in zoom in
  endZoomIn <- PTEN_genomeend + 50*15000 #genome coordinate of last bin included in zoom in
  zoomIn_cnr_nonadj <- which(dat_cnr_nonadj$startcoord>startZoomIn & dat_cnr_nonadj$startcoord<endZoomIn) #which bins start in the zoom in window
  zoomIn_cns_nonadj <- which(dat_cns_nonadj$startcoord>startZoomIn & dat_cns_nonadj$startcoord<endZoomIn) #which segments start in the zoom in window
  segmStartXcoordZoomIn_nonadj <- c(startZoomIn,dat_cns_nonadj$startcoord[zoomIn_cns_nonadj]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn_nonadj <- c(dat_cns_nonadj$endcoord[zoomIn_cns_nonadj-1],endZoomIn) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns_nonadj)>0) {
    segmYcoordZoomIn_nonadj <- dat_cns_nonadj$log2[c((zoomIn_cns_nonadj[1]-1),zoomIn_cns_nonadj)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn_nonadj <- dat_cns_nonadj$log2[which(dat_cns_nonadj$startcoord<=startZoomIn & dat_cns_nonadj$endcoord>=startZoomIn)]
  }
  zoomIn_cnr_adj <- which(dat_cnr_adj$startcoord>startZoomIn & dat_cnr_adj$startcoord<endZoomIn) #which bins start in the zoom in window
  zoomIn_cns_adj <- which(dat_cns_adj$startcoord>startZoomIn & dat_cns_adj$startcoord<endZoomIn) #which segments start in the zoom in window
  segmStartXcoordZoomIn_adj <- c(startZoomIn,dat_cns_adj$startcoord[zoomIn_cns_adj]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn_adj <- c(dat_cns_adj$endcoord[zoomIn_cns_adj-1],endZoomIn) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns_adj)>0) {
    segmYcoordZoomIn_adj <- dat_cns_adj$log2[c((zoomIn_cns_adj[1]-1),zoomIn_cns_adj)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn_adj <- dat_cns_adj$log2[which(dat_cns_adj$startcoord<=startZoomIn & dat_cns_adj$endcoord>=startZoomIn)]
  }
  
  #plot genome wide comparison between pca-adjusted and nonadjusted panel data with flat reference
  title1 <- "panel data, only GC corrected"
  title2 <- "panel data, GC corrected and through pca-clustering"
  title_zoomIn <- "zoom in on PTEN"
  filename <- paste("Bv6_genome_PTEN_compAdjNonadj_",samples[i],".jpeg",sep="")
  xticks <- as.numeric(chr_coord$start)+(as.numeric(chr_coord$end)-as.numeric(chr_coord$start))/2
  
  #plot for tumors with ref flat 
  jpeg(filename,width=9,height=6,quality=90,units="in",res=600)
  layout(matrix(c(1,1,1,2,3,3,3,4), 2, 4, byrow = TRUE))
  par(mar=c(5,4.3,4,2)+0.1)
  #plot nonadjusted bins and segments
  plot(dat_cnr_nonadj$startcoord, dat_cnr_nonadj$log2, xaxt="n", pch=16, col='#00000090', cex=0.3, ylim=c((-3.5),3.5), main=title1,xlab="Chromosome #", ylab="log2 CN ratio",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  segments(dat_cns_nonadj$startcoord, dat_cns_nonadj$log2, dat_cns_nonadj$endcoord, dat_cns_nonadj$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"),cex.axis=1.5)
  abline(v=PTEN_genomestart,lty=2)
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1,cex=1.2)
  #zoom in on PTEN 
  plot(dat_cnr_nonadj$startcoord[zoomIn_cnr_nonadj], dat_cnr_nonadj$log2[zoomIn_cnr_nonadj], xaxt="n", pch=16, col='#000000', cex=0.5, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10 (Mb)", ylab="log2 CN ratio",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  segments(segmStartXcoordZoomIn_nonadj, segmYcoordZoomIn_nonadj, segmEndXcoordZoomIn_nonadj, segmYcoordZoomIn_nonadj, lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7),labels=seq(89,90.5,0.5),cex.axis=1.5)
  text(PTEN_genomeend,3,"PTEN",pos=4,offset=0.2,cex=1.2)
  #plot pca-adjusted bins and segments
  plot(dat_cnr_adj$startcoord, dat_cnr_adj$log2, xaxt="n", pch=16, col='#00000090', cex=0.3, ylim=c((-3.5),3.5), main=title2,xlab="Chromosome #", ylab="log2 CN ratio",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  segments(dat_cns_adj$startcoord, dat_cns_adj$log2, dat_cns_adj$endcoord, dat_cns_adj$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"),cex.axis=1.5)
  abline(v=PTEN_genomestart,lty=2)
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1,cex=1.2)
  #zoom in on PTEN 
  plot(dat_cnr_adj$startcoord[zoomIn_cnr_adj], dat_cnr_adj$log2[zoomIn_cnr_adj], xaxt="n", pch=16, col='#000000', cex=0.5, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10 (Mb)", ylab="log2 CN ratio",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  segments(segmStartXcoordZoomIn_adj, segmYcoordZoomIn_adj, segmEndXcoordZoomIn_adj, segmYcoordZoomIn_adj, lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7),labels=seq(89,90.5,0.5),cex.axis=1.5)
  text(PTEN_genomeend,3,"PTEN",pos=4,offset=0.2,cex=1.2)
  dev.off()  
  
}

write.table(mrpb_nonadj,file=paste(path_CNVkit,"MAPDcomp/mrpb_refFlat_nonadj.txt",sep=""),quote=FALSE,sep="\t")
write.table(mrpb_adj,file=paste(path_CNVkit,"MAPDcomp/mrpb_refFlat_adj.txt",sep=""),quote=FALSE,sep="\t")
write.table(MAPD_nonadj,file=paste(path_CNVkit,"MAPDcomp/MAPD_refFlat_nonadj.txt",sep=""),quote=FALSE,sep="\t")
write.table(MAPD_adj,file=paste(path_CNVkit,"MAPDcomp/MAPD_refFlat_adj.txt",sep=""),quote=FALSE,sep="\t")

