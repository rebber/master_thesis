#Plot the segments and bins genome wide from QDNAseq for wgs data and CNVkit for panel data  with tumor reference for MSI-P2 samples
#Created by Rebecka Bergström on 2016-03-08
library(tools)
path_CNVkit_turef <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT2/100kbATbins/full_tumors/"
path_MSI <- "/home/Crisp/clinseq/MSI-PILOT2/"
setwd(paste(path_CNVkit_turef,"Rplots/segmentPlotsComp/",sep=""))
samples <- unique(basename(file_path_sans_ext(dir(pattern="T-TD3-CB1-panel.bam$",path_MSI,recursive=TRUE)))) #the MSI-P2-samples
samples <- gsub("-TD3-CB1-panel","",samples)

#the files for the panel data
files_bins_turef <- dir(pattern="txt", paste(path_CNVkit_turef,"cnrRcComb",sep=""), recursive = TRUE)
files_segs_turef <- dir(pattern="cns", paste(path_CNVkit_turef,"batchResults",sep=""), recursive = TRUE)
files_rc_targets_turef <- dir(pattern="readcount.targets.txt", paste(path_CNVkit_turef,"readCounts",sep=""),recursive=TRUE)
files_rc_antitargets_turef <- dir(pattern="readcount.antitargets.txt", paste(path_CNVkit_turef,"readCounts",sep=""),recursive=TRUE)

#the files for the wgs data
files_wgs <- dir(pattern="qdnaseq.segments.txt$", path_MSI, recursive = TRUE)

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

#preallocate for MAPD plot
MAPD_total <- matrix(c(samples,rep(NA,2*length(samples))),length(samples),3)
colnames(MAPD_total) <- c("sample","panelTuRef","wgs")
MAPD_targets <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(MAPD_targets) <- c("sample","panelTuRef")
MAPD_antitargets <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(MAPD_antitargets) <- c("sample","panelTuRef")
mrpb_total <- matrix(c(samples,rep(NA,2*length(samples))),length(samples),3)
colnames(mrpb_total) <- c("sample","panelTuRef","wgs")
mrpb_targets <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(mrpb_targets) <- c("sample","panelTuRef")
mrpb_antitargets <- matrix(c(samples,rep(NA,1*length(samples))),length(samples),2)
colnames(mrpb_antitargets) <- c("sample","panelTuRef")

#run for all samples
ptm<-proc.time()
for (i in 1:length(files_bins_turef)) {
  #the panel data with the tumor reference
  dat_bins_turef <- read.table(paste(path_CNVkit_turef,"cnrRcComb/", files_bins_turef[i], sep=""), header = TRUE, as.is=TRUE)
  dat_bins_turef$startcoord <- rep(NA, length(dat_bins_turef$chromosome)) 
  dat_bins_turef$endcoord <- rep(NA, length(dat_bins_turef$chromosome)) 
  dat_segs_turef <- read.table(paste(path_CNVkit_turef,"batchResults/", files_segs_turef[i], sep=""), header = TRUE, as.is=TRUE)
  dat_segs_turef$startcoord <- rep(NA, length(dat_segs_turef$chromosome)) 
  dat_segs_turef$endcoord <- rep(NA, length(dat_segs_turef$chromosome)) 
  #add the genome wide coordinates
  for (j in 1:24) {
    index_bins_turef <- which(dat_bins_turef$chromosome==chr_coord$chromosome[j]) 
    dat_bins_turef$startcoord[index_bins_turef] <- chr_coord$start[j] + dat_bins_turef$start[index_bins_turef] - 1 #add the "genome coordinates" for each bin
    dat_bins_turef$endcoord[index_bins_turef] <- chr_coord$start[j] + dat_bins_turef$end[index_bins_turef] -1 #add the "genome coordinates" for each bin
    index_segs_turef <- which(dat_segs_turef$chromosome==chr_coord$chromosome[j])
    dat_segs_turef$startcoord[index_segs_turef] <- chr_coord$start[j] + dat_segs_turef$start[index_segs_turef] - 1 #add the "genome coordinates" for each segment
    dat_segs_turef$endcoord[index_segs_turef] <- chr_coord$start[j] + dat_segs_turef$end[index_segs_turef] -1 #add the "genome coordinates" for each segment
  }
  
  #the wgs data
  dat_wgs <- read.table(paste(path_MSI, files_wgs[i], sep=""), header = TRUE, as.is=TRUE)
  dat_wgs$startcoord <- rep(NA, length(dat_wgs$chromosome))
  dat_wgs$endcoord <- rep(NA, length(dat_wgs$chromosome))
  #add the genome wide coordinates
  for (j in 1:24) {
    index_wgs <- which(dat_wgs$chromosome==chr_coord$chromosome[j]) 
    dat_wgs$startcoord[index_wgs] <- chr_coord$start[j] + dat_wgs$start[index_wgs] - 1 #add the "genome coordinates" for each bin
    dat_wgs$endcoord[index_wgs] <- chr_coord$start[j] + dat_wgs$end[index_wgs] -1 #add the "genome coordinates" for each bin
  }
  
  #find median read depths per bin (mrpb) and store
  mrpb_total[i,2] <- median(dat_bins_turef$reads)
  mrpb_total[i,3] <- median(dat_wgs$readcount)
  mrpb_targets[i,2] <- median(dat_bins_turef$reads[which(dat_bins_turef$targetType=="T")])
  mrpb_antitargets[i,2] <- median(dat_bins_turef$reads[which(dat_bins_turef$targetType=="AT")])
  print(paste("The median read depth per bin in sample", samples[i], "for panel data with tumor reference is", mrpb_targets[i,2], "reads per bin in targets,", mrpb_antitargets[i,2], "reads per bin in antitargets and", mrpb_total[i,2], "reads per bin in total."))
  print(paste("The median read depth per bin in sample", samples[i], "for wgs data is", mrpb_total[i,3], "reads per bin in total."))
  
  #find the median absolute pairwise difference (MAPD) for the log2CN-ratios and store
  MAPD_total[i,2] <- median(abs(diff(dat_bins_turef$log2)))
  MAPD_total[i,3] <- median(abs(diff(log2(dat_wgs$copynumber))),na.rm=TRUE)
  MAPD_targets[i,2] <- median(abs(diff(dat_bins_turef$log2[which(dat_bins_turef$targetType=="T")])))
  MAPD_antitargets[i,2] <- median(abs(diff(dat_bins_turef$log2[which(dat_bins_turef$targetType=="AT")])))
  
  #the indeces and coordinates for zoom in on PTEN 
  #for wgs
  startZoomIn <- floor(PTEN_genomestart/15000)-50 #index of the first bin included in zoom in
  endZoomIn <- floor(PTEN_genomeend/15000)+50 #index of last bin included in zoom in
  #for the panel
  zoomIn_cnr <- which(dat_bins_turef$startcoord>dat_wgs$startcoord[startZoomIn] & dat_bins_turef$startcoord<dat_wgs$startcoord[endZoomIn]) #which bins start in the zoom in window
  zoomIn_cns <- which(dat_segs_turef$start>dat_wgs$start[startZoomIn] & dat_segs_turef$start<dat_wgs$start[endZoomIn]) #which segments start in the zoom in window
  segmStartXcoordZoomIn <- c(dat_wgs$startcoord[startZoomIn],dat_segs_turef$startcoord[zoomIn_cns]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn <- c(dat_segs_turef$endcoord[zoomIn_cns-1],dat_wgs$endcoord[endZoomIn]) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns)>0) {
    segmYcoordZoomIn <- dat_segs_turef$log2[c((zoomIn_cns[1]-1),zoomIn_cns)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn <- dat_segs_turef$log2[which(dat_segs_turef$startcoord<=dat_wgs$startcoord[startZoomIn] & dat_segs_turef$endcoord>=dat_wgs$startcoord[startZoomIn])]
  }
  
  #plot genome wide comparison between wgs data and panel data with tumor reference
  title1 <- paste("log2 copy number QDNAseq MSI2 ",basename(file_path_sans_ext(files_wgs[i])),sep="")
  title2 <- paste("log2 copy number CNVkit tumor ref MSI2 ",basename(file_path_sans_ext(files_bins_turef[i])),sep="")
  title_wgs_zoomIn <- "zoom in on PTEN, wgs"
  title_panel_zoomIn <- "zoom in on PTEN, panel"
  
  filename <- paste("MSI2_compWgsPanel_genome_PTEN",basename(file_path_sans_ext(gsub("-TD3-CB1-panel","",files_bins_turef[i]))),".pdf",sep="")
  xticks <- as.numeric(chr_coord$start)+(as.numeric(chr_coord$end)-as.numeric(chr_coord$start))/2
  
  pdf(filename,width=20,height=9.8)
  layout(matrix(c(1,1,2,3,3,4), 2, 3, byrow = TRUE)) #plot for wgs and panel in the same graph
  
  #plot wgs
  plot(dat_wgs$startcoord, log2(dat_wgs$copynumber), xaxt="n", pch=16, col='#00000010', cex=0.3, ylim=c((-3.5),3.5), main=title1,xlab="Chromosome #", ylab="log2 CN ratio")
  points(dat_wgs$startcoord, log2(dat_wgs$segmented), pch=16, col="#CD8500",cex=0.6)
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"))
  #zoom in on PTEN
  plot(dat_wgs$startcoord[(startZoomIn):(endZoomIn)],log2(dat_wgs$copynumber[(startZoomIn):(endZoomIn)]), xaxt="n",pch=16, col='#000000', cex=0.3, xlim=(c(dat_wgs$startcoord[(startZoomIn)]-1e4,dat_wgs$startcoord[(endZoomIn)]+1e4)), ylim=c((-3.5),3.5), main=title_wgs_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN wgs")
  segments(dat_wgs$startcoord[(startZoomIn):(endZoomIn)],log2(dat_wgs$segmented[(startZoomIn):(endZoomIn)]),dat_wgs$endcoord[(startZoomIn):(endZoomIn)],log2(dat_wgs$segmented[(startZoomIn):(endZoomIn)]), lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7))
  
  #plot panel 
  plot(dat_bins_turef$startcoord, dat_bins_turef$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title2,xlab="Chromosome #", ylab="log2 CN ratio")
  segments(dat_segs_turef$startcoord, dat_segs_turef$log2, dat_segs_turef$endcoord, dat_segs_turef$log2, lwd=4, col="#CD8500")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=xticks, label=c(1:22,"X","Y"))
  #zoom in on PTEN 
  plot(dat_bins_turef$startcoord[zoomIn_cnr], dat_bins_turef$log2[zoomIn_cnr], xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(dat_wgs$startcoord[(startZoomIn)]-1e4,dat_wgs$startcoord[(endZoomIn)]+1e4)), ylim=c((-3.5),3.5), main=title_panel_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN panel")
  segments(segmStartXcoordZoomIn,segmYcoordZoomIn,segmEndXcoordZoomIn,segmYcoordZoomIn, lwd=4, col="#CD8500")
  abline(v=c(PTEN_genomestart,PTEN_genomeend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=as.numeric(chr_coord$end[9])+seq(8.9e7,9.05e7,0.05e7) ,labels=seq(8.9e7,9.05e7,0.05e7))
  
  dev.off()
  
}

proc.time()-ptm

#plot MAPD plots
setwd(paste(path_CNVkit_turef,"Rplots/",sep=""))
reads=1:2000
mapd<-c()
for (i in 1:2000) mapd[i]=median(abs(diff(log2(rpois(10000,reads[i])/reads[i]))))
filename_MAPD <- "MAPD_MSI2.pdf"
pdf(filename_MAPD)
layout(matrix(1,1,1,byrow=TRUE))
plot(reads,mapd,xlab='median #reads per bin',ylab='median absolute pairwise difference (log2CNratio)',type='l',main='Minimum stochastic noise',log="x",ylim=c(0,0.8))
points(reads, mapd*sqrt(2),type="l",lty=2)
#points(as.double(mrpb_total[,2]),MAPD_total[,2],col="red")
points(as.double(mrpb_total[,3]),MAPD_total[,3],col="blue")
points(as.double(mrpb_targets[,2]),MAPD_targets[,2],col="darkorchid4")
points(as.double(mrpb_antitargets[,2]),MAPD_antitargets[,2],col="darkgreen")
legend("topright", c("all bins, wgs","targets, with tumor ref","antitargets, with tumor ref","Poisson distribution","Poisson*sqrt(2)"),col=c("blue","darkorchid4","darkgreen","black","black"),pch=c(rep(1,3),NA,NA),lty=c(rep(NA,3),1,2))
dev.off()

filename_MAPDvsMAPD <- "MAPDvsMAPD_mrpbVSmrpb.pdf"
pdf(filename_MAPDvsMAPD)
plot(MAPD_targets[,2],MAPD_antitargets[,2])
plot(mrpb_targets[,2],mrpb_antitargets[,2])
dev.off()

filename_density <- "densities_MSI2.pdf"
pdf(filename_density)
plot(density(dat_bins_turef$reads)) #,xlim=c(-100,2000)
plot(density(dat_bins_turef$end-dat_bins_turef$start))
dev.off()

quantile(dat_bins_turef$end-dat_bins_turef$start,probs=seq(0,1,0.1))
shortBins <- which((dat_bins_turef$end-dat_bins_turef$start)<50000)
median((dat_bins_turef$end[shortBins]-dat_bins_turef$start[shortBins]))
quantile((dat_bins_turef$end[shortBins]-dat_bins_turef$start[shortBins]),probs=seq(0,1,0.1))
AT_bins_ind <- which(dat_bins_turef$targetType=="AT")
mean(dat_bins_turef$end[AT_bins_ind]-dat_bins_turef$start[AT_bins_ind])
median(dat_bins_turef$end[AT_bins_ind]-dat_bins_turef$start[AT_bins_ind])
length(which((dat_bins_turef$end[AT_bins_ind]-dat_bins_turef$start[AT_bins_ind])>1e6))

View(cbind(mrpb_total,mrpb_targets[,2],mrpb_antitargets[,2],signif(as.numeric(mrpb_targets[,2])/as.numeric(mrpb_antitargets[,2]),digits=3)))
median(as.numeric(mrpb_targets[,2])/as.numeric(mrpb_antitargets[,2]))
mean(as.numeric(mrpb_targets[,2])/as.numeric(mrpb_antitargets[,2]))
median(as.numeric(mrpb_total[,2]))
median(as.numeric(mrpb_total[,3]))
median(as.numeric(mrpb_targets[,2]))
median(as.numeric(mrpb_antitargets[,2]))

#layout(1)
#qqnorm(dat_cnr_blref$log2[which(dat_cnr_blref$chromosome!="Y")])
#qqline(dat_cnr_blref$log2)
#qqnorm(dat_cnr_turef$log2)
#qqline(dat_cnr_turef$log2)
#qqplot(dat_cnr_blref$log2,dat_cnr_turef$log2,ylim=c(-3.5,3.5),xlim=c(-3.5,3.5))
#abline(0,1)

#qqplot(dat_cnr_blref$log2[which(dat_cnr_blref$chromosome!="Y")],dat_cns_blref$log2[which(dat_cnr_blref$chromosome!="Y")])
#abline(0,1)



