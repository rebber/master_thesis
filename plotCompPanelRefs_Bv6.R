#Plot the segments and bins from CNVkit for panel data with blood reference vs with tumor reference for Bv6 testsamples 
library(tools)
path_CNVkit_blref <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/ref5norm/"
path_CNVkit_turef <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/ref10tum/"
setwd(paste(path_CNVkit_turef,"Rplots/refComp/",sep=""))
load("~/cnvkit/BREASTv6_all/test50/test50") #load test set samples
samples <- testset[,1]

files_cnr_blref <- dir(pattern="cnr", paste(path_CNVkit_blref,"batchResults",sep=""), recursive = TRUE)
files_cns_blref <- dir(pattern="cns", paste(path_CNVkit_blref,"batchResults",sep=""), recursive = TRUE)
files_cnr_turef <- dir(pattern="cnr", paste(path_CNVkit_turef,"batchResults",sep=""), recursive = TRUE)
files_cns_turef <- dir(pattern="cns", paste(path_CNVkit_turef,"batchResults",sep=""), recursive = TRUE)

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


#run for all samples
ptm<-proc.time()
for (i in 1:length(files_cnr_blref)) {
  #the data from the blood reference
  dat_cnr_blref <- read.table(paste(path_CNVkit_blref,"batchResults/", files_cnr_blref[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_blref$startcoord <- rep(NA, length(dat_cnr_blref$chromosome)) 
  dat_cnr_blref$endcoord <- rep(NA, length(dat_cnr_blref$chromosome)) 
  dat_cns_blref <- read.table(paste(path_CNVkit_blref,"batchResults/", files_cns_blref[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cns_blref$startcoord <- rep(NA, length(dat_cns_blref$chromosome)) 
  dat_cns_blref$endcoord <- rep(NA, length(dat_cns_blref$chromosome))  
  #add the genome wide coordinates for each bin
  for (j in 1:24) {
    index_cnr_blref <- which(dat_cnr_blref$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_blref$startcoord[index_cnr_blref] <- chr_coord$start[j] + dat_cnr_blref$start[index_cnr_blref] - 1 #add the "genome coordinates" for each bin
    dat_cnr_blref$endcoord[index_cnr_blref] <- chr_coord$start[j] + dat_cnr_blref$end[index_cnr_blref] -1 #add the "genome coordinates" for each bin
    index_cns_blref <- which(dat_cns_blref$chromosome==chr_coord$chromosome[j])
    dat_cns_blref$startcoord[index_cns_blref] <- chr_coord$start[j] + dat_cns_blref$start[index_cns_blref] - 1 #add the "genome coordinates" for each segment
    dat_cns_blref$endcoord[index_cns_blref] <- chr_coord$start[j] + dat_cns_blref$end[index_cns_blref] -1 #add the "genome coordinates" for each segment
  }
  
  #the data from the tumor reference
  dat_cnr_turef <- read.table(paste(path_CNVkit_turef,"batchResults/", files_cnr_turef[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr_turef$startcoord <- rep(NA, length(dat_cnr_turef$chromosome)) 
  dat_cnr_turef$endcoord <- rep(NA, length(dat_cnr_turef$chromosome)) 
  dat_cns_turef <- read.table(paste(path_CNVkit_turef,"batchResults/", files_cns_turef[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cns_turef$startcoord <- rep(NA, length(dat_cns_turef$chromosome)) 
  dat_cns_turef$endcoord <- rep(NA, length(dat_cns_turef$chromosome)) 
  #add the genome wide coordinates
  for (j in 1:24) {
    index_cnr_turef <- which(dat_cnr_turef$chromosome==chr_coord$chromosome[j]) 
    dat_cnr_turef$startcoord[index_cnr_turef] <- chr_coord$start[j] + dat_cnr_turef$start[index_cnr_turef] - 1 #add the "genome coordinates" for each bin
    dat_cnr_turef$endcoord[index_cnr_turef] <- chr_coord$start[j] + dat_cnr_turef$end[index_cnr_turef] -1 #add the "genome coordinates" for each bin
    index_cns_turef <- which(dat_cns_turef$chromosome==chr_coord$chromosome[j])
    dat_cns_turef$startcoord[index_cns_turef] <- chr_coord$start[j] + dat_cns_turef$start[index_cns_turef] - 1 #add the "genome coordinates" for each segment
    dat_cns_turef$endcoord[index_cns_turef] <- chr_coord$start[j] + dat_cns_turef$end[index_cns_turef] -1 #add the "genome coordinates" for each segment
  }
  
  #the indeces and coordinates for zoom in on PTEN 
  startZoomIn <- PTEN_genomestart - 50*15000 #genome coordinate of the first bin included in zoom in
  endZoomIn <- PTEN_genomeend + 50*15000 #genome coordinate of last bin included in zoom in
  zoomIn_cnr_blref <- which(dat_cnr_blref$startcoord>startZoomIn & dat_cnr_blref$startcoord<endZoomIn) #which bins start in the zoom in window
  zoomIn_cns_blref <- which(dat_cns_blref$startcoord>startZoomIn & dat_cns_blref$startcoord<endZoomIn) #which segments start in the zoom in window
  segmStartXcoordZoomIn_blref <- c(startZoomIn,dat_cns_blref$startcoord[zoomIn_cns_blref]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn_blref <- c(dat_cns_blref$endcoord[zoomIn_cns_blref-1],endZoomIn) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns_blref)>0) {
    segmYcoordZoomIn_blref <- dat_cns_blref$log2[c((zoomIn_cns_blref[1]-1),zoomIn_cns_blref)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn_blref <- dat_cns_blref$log2[which(dat_cns_blref$startcoord<=startZoomIn & dat_cns_blref$endcoord>=startZoomIn)]
  }
  zoomIn_cnr_turef <- which(dat_cnr_turef$startcoord>startZoomIn & dat_cnr_turef$startcoord<endZoomIn) #which bins start in the zoom in window
  zoomIn_cns_turef <- which(dat_cns_turef$startcoord>startZoomIn & dat_cns_turef$startcoord<endZoomIn) #which segments start in the zoom in window
  segmStartXcoordZoomIn_turef <- c(startZoomIn,dat_cns_turef$startcoord[zoomIn_cns_turef]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn_turef <- c(dat_cns_turef$endcoord[zoomIn_cns_turef-1],endZoomIn) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns_turef)>0) {
    segmYcoordZoomIn_turef <- dat_cns_turef$log2[c((zoomIn_cns_turef[1]-1),zoomIn_cns_turef)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn_turef <- dat_cns_turef$log2[which(dat_cns_turef$startcoord<=startZoomIn & dat_cns_turef$endcoord>=startZoomIn)]
  }
  
  #plot genome wide comparison between panel data with blood reference and tumor reference
  title1 <- "panel data with 5 normals as reference"
  title2 <- "panel data with 10 tumors as reference"
  title_zoomIn <- "zoom in on PTEN"
  filename <- paste("Bv6_compPanelRefs_",samples[i],".jpeg",sep="")
  xticks <- as.numeric(chr_coord$start)+(as.numeric(chr_coord$end)-as.numeric(chr_coord$start))/2
  
  jpeg(filename,width=9,height=6,quality=90,units="in",res=600)
  layout(matrix(c(1,1,1,2,3,3,3,4), 2, 4, byrow = TRUE))
  par(mar=c(5,4.3,4,2)+0.1)
  plot(dat_cnr_blref$startcoord, dat_cnr_blref$log2, xaxt="n", pch=16, col='#00000090', cex=0.3, ylim=c((-3.5),3.5), main=title1,xlab="Chromosome #", ylab="log2 CN ratio",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  segments(dat_cns_blref$startcoord, dat_cns_blref$log2, dat_cns_blref$endcoord, dat_cns_blref$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"),cex.axis=1.5)
  abline(v=PTEN_genomestart,lty=2)
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1,cex=1.2)
  #zoom in on PTEN 
  plot(dat_cnr_blref$startcoord[zoomIn_cnr_blref], dat_cnr_blref$log2[zoomIn_cnr_blref], xaxt="n", pch=16, col='#000000', cex=0.5, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10 (Mb)", ylab="log2 CN ratio",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  segments(segmStartXcoordZoomIn_blref, segmYcoordZoomIn_blref, segmEndXcoordZoomIn_blref, segmYcoordZoomIn_blref, lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7),labels=seq(89,90.5,0.5),cex.axis=1.5)
  text(PTEN_genomeend,3,"PTEN",pos=4,offset=0.2,cex=1.2)
  
  plot(dat_cnr_turef$startcoord, dat_cnr_turef$log2, xaxt="n", pch=16, col='#00000090', cex=0.3, ylim=c((-3.5),3.5), main=title2,xlab="Chromosome #", ylab="log2 CN ratio",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  segments(dat_cns_turef$startcoord, dat_cns_turef$log2, dat_cns_turef$endcoord, dat_cns_turef$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"),cex.axis=1.5)
  abline(v=PTEN_genomestart,lty=2)
  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1,cex=1.2)
  #zoom in on PTEN 
  plot(dat_cnr_turef$startcoord[zoomIn_cnr_turef], dat_cnr_turef$log2[zoomIn_cnr_turef], xaxt="n", pch=16, col='#000000', cex=0.5, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10 (Mb)", ylab="log2 CN ratio",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
  segments(segmStartXcoordZoomIn_turef, segmYcoordZoomIn_turef, segmEndXcoordZoomIn_turef, segmYcoordZoomIn_turef, lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7),labels=seq(89,90.5,0.5),cex.axis=1.5)
  text(PTEN_genomeend,3,"PTEN",pos=4,offset=0.2,cex=1.2)
  
  dev.off()
  
}

proc.time()-ptm