# Plot segments and bins (wgs & panel) for one sample, where CNVkit segmentation had less stringent cut-off
#Created by Rebecka Bergstrom at 2016-03-11

library(tools)
path_CNVkit_turef <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT2/100kbATbins/full_tumors/"
path_MSI <- "/home/Crisp/clinseq/MSI-PILOT2/"
pname <- "p1e-2" #the p-value cut-off
sample <- "ex1" #the sample I wanna look at (exchange ex1 for actual sample ID)
setwd(paste(path_CNVkit_turef,pname,"/Rplots/segmentPlotsComp/",sep=""))
samples <- unique(basename(file_path_sans_ext(dir(pattern="T-TD3-CB1-panel.bam$",path_MSI,recursive=TRUE)))) #the MSI-P2-samples
samples <- gsub("-TD3-CB1-panel","",samples)


#the files for the panel data
files_bins_turef <- dir(pattern="txt", paste(path_CNVkit_turef,"cnrRcComb",sep=""), recursive = TRUE)
files_segs_turef <- dir(pattern="cns", paste(path_CNVkit_turef,pname,"/cns/",sep=""), recursive = TRUE)
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


i<- grep(sample, files_bins_turef)

#the panel data with the tumor reference
dat_bins_turef <- read.table(paste(path_CNVkit_turef,"cnrRcComb/", files_bins_turef[i], sep=""), header = TRUE, as.is=TRUE)
dat_bins_turef$startcoord <- rep(NA, length(dat_bins_turef$chromosome)) 
dat_bins_turef$endcoord <- rep(NA, length(dat_bins_turef$chromosome)) 
dat_segs_turef <- read.table(paste(path_CNVkit_turef,pname,"/cns/", files_segs_turef[i], sep=""), header = TRUE, as.is=TRUE)
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
title2 <- paste("log2 copy number CNVkit tumor ref MSI2, p-value segments ",pname,basename(file_path_sans_ext(files_bins_turef[i])),sep="")
title_wgs_zoomIn <- "zoom in on PTEN, wgs"
title_panel_zoomIn <- "zoom in on PTEN, panel"

filename <- paste("MSI2_compWgsPanel_genome_PTEN_",pname,"_",basename(file_path_sans_ext(gsub("-TD3-CB1-panel","",files_bins_turef[i]))),".pdf",sep="")
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


#Test differences between segment for p=1e-4 and p=1e-3
p4 <- read.table("/home/rebecka.bergstrom/cnvkit/MSI-PILOT2/100kbATbins/full_tumors/batchResults/ex1T-TD3-CB1-panel.cns",header=TRUE,as.is=TRUE)
p_higher <- read.table(paste("/home/rebecka.bergstrom/cnvkit/MSI-PILOT2/100kbATbins/full_tumors/",pname,"/cns/ex1T-TD3-CB1-panel_",pname,".cns",sep=""),header=TRUE,as.is=TRUE)

#The comments are valid for p=1e-3
identical(p4,p_higher) #They are not completely identical
identical(p4$start,p_higher$start) #But identical start-
identical(p4$end,p_higher$end) #and end-positions
identical(p4$log2,p_higher$log2) #and log2-CN values
identical(p4$gene,p_higher$gene) #and genes
identical(p4$probes,p_higher$probes) #and probes
identical(p4$chromosome,p_higher$chromosome) #and chromosome #
identical(p4$weight,p_higher$weight) #However the weights are not identical
median(p4$weight-p_higher$weight)
quantile(p4$weight-p_higher$weight)
which(p4$weight!=p_higher$weight) #only 2 of 70 segments (seg 13 and 70) differ between the to cases

#For p=1e-2 there are 3 more segments than for p=1e-4

