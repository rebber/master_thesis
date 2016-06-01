# Compare segmentation and bins for Bv6 test samples with ref of 5 normals
#Created by Rebecka Bergström on 2016-04-20

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/ref5norm/"
path_wgs <- "/home/Crisp/clinseq/from_rasta/BREASTv6/"
load("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/test50/test50")
samples <- testset[,1]
#files_cnr <- paste("batchResults/",samples,".cnr",sep="")
files_cnr <- paste("cnrRcComb/addedSmooth/",samples,".txt",sep="")
files_cns <- paste("batchResults/",samples,".cns",sep="")
files_wgs <- rep(NA, length(samples))
for (i in 1:length(samples)) {
  files_wgs[i] <- dir(pattern=paste(gsub("_panel_v1","",samples[i]),".*qdnaseq.txt$",sep=""), path=path_wgs,recursive=TRUE)
}

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


#read files and plot genomewide for all samples
for (i in 1:length(samples)) {
  dat_cnr <- read.table(paste(path_CNVkit, files_cnr[i],sep=""),header=TRUE,as.is=TRUE)
  dat_cnr$startcoord <- rep(NA,length(dat_cnr$chromosome))
  dat_cnr$endcoord <- rep(NA,length(dat_cnr$chromosome))
  dat_cns <- read.table(paste(path_CNVkit, files_cns[i],sep=""),header=TRUE,as.is=TRUE)  
  dat_cns$startcoord <- rep(NA,length(dat_cns$chromosome))
  dat_cns$endcoord <- rep(NA,length(dat_cns$chromosome))
  dat_wgs <- read.table(paste(path_wgs,files_wgs[i],sep=""),header=TRUE,as.is=TRUE)
  dat_wgs$startcoord <- rep(NA,length(dat_wgs$chromosome))
  dat_wgs$endcoord <- rep(NA,length(dat_wgs$chromosome))
  
  #add the genome wide coordinates
  for (j in 1:24) {
    #add the "genome coordinates" for each bin in panel 
    index_cnr <- which(dat_cnr$chromosome==chr_coord$chromosome[j]) 
    dat_cnr$startcoord[index_cnr] <- chr_coord$start[j] + dat_cnr$start[index_cnr] - 1 
    dat_cnr$endcoord[index_cnr] <- chr_coord$start[j] + dat_cnr$end[index_cnr] -1 
    #add the "genome coordinates" for each segment in panel
    index_cns <- which(dat_cns$chromosome==chr_coord$chromosome[j])
    dat_cns$startcoord[index_cns] <- chr_coord$start[j] + dat_cns$start[index_cns] - 1 
    dat_cns$endcoord[index_cns] <- chr_coord$start[j] + dat_cns$end[index_cns] -1 
    #add the "genome coordinates" for each bin in the wgs data
    index_wgs <- which(dat_wgs$chromosome==chr_coord$chromosome[j])
    dat_wgs$startcoord[index_wgs] <- chr_coord$start[j] + dat_wgs$start[index_wgs] - 1 
    dat_wgs$endcoord[index_wgs] <- chr_coord$start[j] + dat_wgs$end[index_wgs] -1 
  }
  
  
  #the indeces and coordinates for zoom in on PTEN 
  startZoomIn <- PTEN_genomestart - 50*15000 #genome coordinate of the first bin included in zoom in
  endZoomIn <- PTEN_genomeend + 50*15000 #genome coordinate of last bin included in zoom in
  zoomIn_cnr <- which(dat_cnr$startcoord>=startZoomIn & dat_cnr$startcoord<endZoomIn) #which bins start in the zoom in window
  zoomIn_cns <- which(dat_cns$startcoord>=startZoomIn & dat_cns$startcoord<endZoomIn) #which segments start in the zoom in window
  segmStartXcoordZoomIn <- c(startZoomIn,dat_cns$startcoord[zoomIn_cns]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn <- c(dat_cns$endcoord[zoomIn_cns-1],endZoomIn) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns)>0) {
    segmYcoordZoomIn <- dat_cns$log2[c((zoomIn_cns[1]-1),zoomIn_cns)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn <- dat_cns$log2[which(dat_cns$startcoord<=startZoomIn & dat_cns$endcoord>=startZoomIn)]
  }
  zoomIn_wgs <- which(dat_wgs$startcoord>=startZoomIn & dat_wgs$startcoord<endZoomIn) #which bins start in the zoom in window
  
  #plot genome wide comparison between wgs (QDNAseq) data and panel (CNVkit) data
  title1 <- "wgs data"
  title2 <- "panel data"
  title_zoomIn <- "zoom in on PTEN"
  filename <- paste("Bv6_genome_PTEN_compWgsPanel_",samples[i],".jpeg",sep="")
  xticks <- as.numeric(chr_coord$start)+(as.numeric(chr_coord$end)-as.numeric(chr_coord$start))/2
  
  #plot genome wide 
  setwd(paste(path_CNVkit,"Rplots/segmentsComp",sep=""))
  jpeg(filename,width=9,height=9,quality=90,units="in",res=600)
  #  layout(matrix(c(1,1,1,2,3,3,3,4), 2, 4, byrow = TRUE)) 
  layout(matrix(c(1,2),nrow=2,ncol=1,byrow=TRUE))
  #plot wgs bins and segments
  plot(dat_wgs$startcoord, log2(dat_wgs$copynumber), xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title1,xlab="Chromosome #", ylab="log2 CN ratio",cex.main=1,cex.axis=1,cex.lab=1)
  segments(dat_wgs$startcoord, log2(dat_wgs$segmented), dat_wgs$endcoord, log2(dat_wgs$segmented), lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"),cex.axis=1)
  #  abline(v=PTEN_genomestart,lty=2)
  #  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1,cex=1.2)
  #zoom in on PTEN 
  #  plot(dat_wgs$startcoord[zoomIn_wgs], log2(dat_wgs$copynumber[zoomIn_wgs]), xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10",cex.main=1,cex.axis=1,cex.lab=1, ylab="log2 CN ratio")
  #  segments(dat_wgs$startcoord[zoomIn_wgs], log2(dat_wgs$segmented[zoomIn_wgs]), dat_wgs$endcoord[zoomIn_wgs], log2(dat_wgs$segmented[zoomIn_wgs]), lwd=4, col="#CD8500")
  # abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  #  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7),cex.axis=1)
  #  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1,cex=1.2)
  #plot panel bins and segments
  plot(dat_cnr$startcoord[which(dat_cnr$targetType=="T")], dat_cnr$log2[which(dat_cnr$targetType=="T")], xaxt="n", pch=16, col='#00000090', cex=0.3, ylim=c((-3.5),3.5), main=title2,xlab="Chromosome #", ylab="log2 CN ratio",cex.main=1,cex.axis=1,cex.lab=1)
  segments(dat_cns$startcoord, dat_cns$log2, dat_cns$endcoord, dat_cns$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"),cex.axis=1)
  #  abline(v=PTEN_genomestart,lty=2)
  #  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1,cex=1.2)
  #zoom in on PTEN 
  #  plot(dat_cnr$startcoord[zoomIn_cnr], dat_cnr$log2[zoomIn_cnr], xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(startZoomIn-1e4,endZoomIn+1e4)), ylim=c((-3.5),3.5), main=title_zoomIn,xlab="Position on Chr 10",cex.main=1,cex.axis=1,cex.lab=1, ylab="log2 CN ratio")
  #  segments(segmStartXcoordZoomIn, segmYcoordZoomIn, segmEndXcoordZoomIn, segmYcoordZoomIn, lwd=4, col="#CD8500")
  #  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  #  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7),cex.axis=1)
  #  text(PTEN_genomestart,3,"PTEN",pos=4,offset=0.1,cex=1.2)
  dev.off() 
  
}
