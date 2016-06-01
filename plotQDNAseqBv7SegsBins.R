#Plot the segments and bins from QDNAseq for wgs data from BREASTv7 samples
setwd("/home/rebecka.bergstrom/BREASTv7/plotsSegmentsBins/")
library(tools)

files <- dir(pattern="T.*qdnaseq.txt", "/home/Crisp/clinseq/BREASTv7/", recursive = TRUE)

for (i in 1:length(files)) {
  dat <- read.table(paste("/home/Crisp/clinseq/BREASTv7/", files[i], sep=""), header = TRUE, as.is=TRUE)
  chr<-dat$chromosome #the chromosome numbers
  chr_ind<-c()
  for (j in 1:22) { #get the indeces where each chromosome starts in chr
    chr_ind[j] <- min(which(chr==j))
  }
  chr_ind <- c(chr_ind, min(which(chr=="X")), min(which(chr=="Y"))) #add indeces for X and Y
  title <- paste("BREASTv7 log2 copy number",basename(file_path_sans_ext(files[i])))
  filename=paste("Bv7_segmentsBins_",basename(file_path_sans_ext(files[i])),".pdf",sep="")
  pdf(filename,width=20,height=9.8)
  plot(log2(dat$copynumber), xaxt="n", pch=16, col='#00000010', cex=0.3, ylim=c((-5),5), main=title,xlab="Chromosome #",cex.main=0.7, ylab="log2 CN")
  points(log2(dat$segmented), cex=0.6, pch=16, col="#CD8500")
  abline(v=c(chr_ind, max(which(chr=="Y"))), h=0, lwd=0.5) #vertical lines for the chromosomes
  axis(1,at=chr_ind, label=c(1:22,"X","Y"))
  dev.off()
}


