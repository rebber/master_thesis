#Plot the segments and bins from CNVkit for panel data for MSI samples
setwd("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/Rplots/segmentPlotsPanel/")

#Get the lengths of all chromosomes and find at which genomic coordinates the chromosomes start and end
files_cnr <- dir(pattern="cnr", "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResults/", recursive = TRUE)
files_cnr <- files_cnr[-grep("broken_pipe",files_cnr)]
files_cns <- dir(pattern="cns", "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResults/", recursive = TRUE)
files_cns <- files_cns[-grep("broken_pipe",files_cns)]
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

for (i in 1:length(files_cnr)) {
  dat_cnr <- read.table(paste("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResults/", files_cnr[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr$startcoord <- rep(NA, length(dat_cnr$chromosome)) 
  dat_cnr$endcoord <- rep(NA, length(dat_cnr$chromosome)) 
  dat_cns <- read.table(paste("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResults/", files_cns[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cns$startcoord <- rep(NA, length(dat_cns$chromosome)) 
  dat_cns$endcoord <- rep(NA, length(dat_cns$chromosome)) 
  for (j in 1:24) {
    index_cnr <- which(dat_cnr$chromosome==chr_coord$chromosome[j]) 
    dat_cnr$startcoord[index_cnr] <- chr_coord$start[j] + dat_cnr$start[index_cnr] - 1 #add the "genome coordinates" for each bin
    dat_cnr$endcoord[index_cnr] <- chr_coord$start[j] + dat_cnr$end[index_cnr] -1 #add the "genome coordinates" for each bin
    index_cns <- which(dat_cns$chromosome==chr_coord$chromosome[j])
    dat_cns$startcoord[index_cns] <- chr_coord$start[j] + dat_cns$start[index_cns] - 1 #add the "genome coordinates" for each segment
    dat_cns$endcoord[index_cns] <- chr_coord$start[j] + dat_cns$end[index_cns] -1 #add the "genome coordinates" for each segment
  }
  
  title <- paste("MSI log2 copy number ",basename(file_path_sans_ext(files_cnr[i])),".cnvkit",sep="")
  filename=paste("MSI_segmentsBins_",basename(file_path_sans_ext(files_cnr[i])),".cnvkit",".pdf",sep="")
  pdf(filename,width=17,height=5)
  plot(dat_cnr$startcoord, dat_cnr$log2, xaxt="n", pch=16, col='#00000030', cex=0.3, ylim=c((-3.5),3.5), main=title,xlab="Chromosome #",cex.main=0.7, ylab="log2 CN")
  segments(dat_cns$startcoord, dat_cns$log2, dat_cns$endcoord, dat_cns$log2, lwd=4, col="#66CD00")
  abline(v=c(chr_coord$start,chr_coord$end[24]), h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  axis(1, at=chr_coord$start, label=c(1:22,"X","Y"))
  dev.off()
}
