#Plot chromosome 10 for wgs QDNAseqdata for samples with PTEN-CN < 0.75
setwd("/home/rebecka.bergstrom/BREASTv7_PTEN_del/QDNAseqWGS")
library(tools)

#Construct PTEN_bed, PTEN_deletion and numBins with script PTEN_CN_BREASTv7.R
files <- PTEN_deletion
PTEN_genestart <- PTEN_bed$start
PTEN_segs <- as.data.frame(cbind(files,matrix(rep(NA,length(files)*numBins),ncol=numBins)),stringsAsFactors=FALSE)

for (i in 1:length(files)) {
  dat <- read.table(paste("/home/Crisp/clinseq/BREASTv7/", files[i], sep=""), header = TRUE, as.is=TRUE)
  chr<-dat$chromosome #the chromosome numbers
  chr_ind<-c()
  chr10_start <- min(which(chr==10))
  chr10_end <- max(which(chr==10))
  title <- paste("BREASTv7 log2 copy number chr 10",basename(file_path_sans_ext(files[i])))
  filename=paste("BREASTv7_chr10_segmentsBins_",basename(file_path_sans_ext(files[i])),".pdf",sep="")
  pdf(filename,width=17,height=5)
  plot(dat$start[chr10_start:chr10_end], log2(dat$copynumber[chr10_start:chr10_end]), xaxt="n", pch=16, col='#00000090', cex=0.3, ylim=c((-3.5),3.5), main=title,xlab="Position on Chr 10",cex.main=0.7, ylab="log2 CN")
  points(dat$start[chr10_start:chr10_end], log2(dat$segmented[chr10_start:chr10_end]), cex=0.6, pch=16, col="#66CD00")
  abline(v=c(PTEN_genestart), h=0, lwd=0.5) #vertical lines for the chromosomes
  axis(1,at=c(dat$end[floor(chr10_start:chr10_end/1000)*1000]))
  dev.off()
  
  PTEN_segs[i,2:10] <- dat$segmented[PTEN_startind:PTEN_endind]
}

PTEN_segs_medians <- cbind(files,rep(NA,length(files)))
for (i in 1:length(PTEN_segs$files)) {
  PTEN_segs_medians[i,2] <- median(as.matrix(PTEN_segs[i,2:10]))
}

quantile(log2(as.numeric(as.matrix(PTEN_segs[,2:10]))),probs=seq(0,1,0.125))
sum(log2(as.numeric(as.matrix(PTEN_segs[,2])))<=-1)
quantile(log2(as.numeric(PTEN_segs_medians[,2])))
sum(log2(as.numeric(PTEN_segs_medians[,2]))<=-0.6)
PTEN_segs_medians[which(log2(as.numeric(PTEN_segs_medians[,2]))<=-0.6),1]
log2(0.66)
2^-1
