#Plot chromosome 17 for wgs QDNAseqdata (BREASTv7) for samples with log2(ERBB2-CN) > 3
setwd("/home/rebecka.bergstrom/BREASTv7_PTEN_del/ERBB2amp/QDNAseqWGS")
library(tools)

#Construct ERBB2_bed, ERBB2_amplification, ERBB2_startind, ERBB2_endtind and numBins_ERBB2 with script ERBB2_CN_BREASTv7.R
files <- ERBB2_amplification
ERBB2_genestart <- ERBB2_bed$start
ERBB2_geneend <- ERBB2_bed$end
ERBB2_segs <- as.data.frame(cbind(files,matrix(rep(NA,length(files)*numBins_ERBB2),ncol=numBins_ERBB2)),stringsAsFactors=FALSE)

for (i in 1:length(files)) {
  dat <- read.table(paste("/home/Crisp/clinseq/BREASTv7/", files[i], sep=""), header = TRUE, as.is=TRUE)
  chr<-dat$chromosome #the chromosome numbers
  chr_ind<-c()
  chr17_start <- min(which(chr==17))
  chr17_end <- max(which(chr==17))
  title <- paste("BREASTv7 log2 copy number chr 17",basename(file_path_sans_ext(files[i])))
  filename=paste("BREASTv7_chr17_segmentsBins_",basename(file_path_sans_ext(files[i])),".pdf",sep="")
  pdf(filename,width=18,height=5)
  layout(matrix(c(1,1,2), 1, 3, byrow = TRUE))
  plot(dat$start[chr17_start:chr17_end], log2(dat$copynumber[chr17_start:chr17_end]), xaxt="n", pch=16, col='#00000090', cex=0.3, ylim=c((-3.5),6), main=title,xlab="Position on Chr 17",cex.main=0.7, ylab="log2 CN")
  points(dat$start[chr17_start:chr17_end], log2(dat$segmented[chr17_start:chr17_end]), cex=0.6, pch=16, col="#66CD00")
  abline(v=c(ERBB2_genestart), h=0, lwd=0.5) #vertical line for ERBB2 and horizontal line for neutral CN
  axis(1,at=c(dat$end[floor(chr17_start:chr17_end/1000)*1000]))
  #also plot zoom in on ERBB2
  plot(dat$start[(ERBB2_startind-100):(ERBB2_endind+100)], log2(dat$copynumber[(ERBB2_startind-100):(ERBB2_endind+100)]), xaxt="n", pch=16, col='#000000', cex=0.3, ylim=c((-3.5),6), main=title,xlab="Position on Chr 17",cex.main=0.7, ylab="log2 CN")
  points(dat$start[(ERBB2_startind-100):(ERBB2_endind+100)], log2(dat$segmented[(ERBB2_startind-100):(ERBB2_endind+100)]), cex=0.6, pch=16, col="#66CD00")
  abline(v=c(ERBB2_genestart,ERBB2_geneend), h=0, lwd=0.5) #vertical line for ERBB2 and horizontal line for neutral CN
  axis(1,at=c(dat$end[floor(((ERBB2_startind-100):(ERBB2_endind+100))/50)*50]))
  dev.off()
  
  ERBB2_segs[i,2:(numBins_ERBB2+1)] <- dat$segmented[ERBB2_startind:ERBB2_endind]
}

ERBB2_segs_medians <- cbind(files,rep(NA,length(files)))
for (i in 1:length(ERBB2_segs$files)) {
  ERBB2_segs_medians[i,2] <- median(as.matrix(ERBB2_segs[i,2:(numBins_ERBB2+1)]))
}

quantile(log2(as.numeric(as.matrix(ERBB2_segs[,2:(numBins_ERBB2+1)]))),probs=seq(0,1,0.125))
sum(log2(as.numeric(as.matrix(ERBB2_segs[,2])))<=-1)
quantile(log2(as.numeric(ERBB2_segs_medians[,2])))
sum(log2(as.numeric(ERBB2_segs_medians[,2]))<=-0.6)
ERBB2_segs_medians[which(log2(as.numeric(ERBB2_segs_medians[,2]))<=-0.6),1]
log2(0.66)
2^-1
