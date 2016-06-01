#Plot CNVkit data
library(tools)

files_cnr <- dir(pattern="cnr", "/home/rebecka.bergstrom/cnvkit/1stRunBreast/batchRun/", recursive = TRUE)
files_cns <- dir(pattern="cns", "/home/rebecka.bergstrom/cnvkit/1stRunBreast/batchRun/", recursive = TRUE)

for (i in 1:length(files_cnr)) {
  cnr <- read.table(paste("/home/rebecka.bergstrom/cnvkit/1stRunBreast/batchRun/", files_cnr[i], sep=""), header = TRUE, as.is=TRUE)
  binlengths<-cnr$end-cnr$start
  quantile(binlengths)
  binmiddles <- (cnr$end-cnr$start)/2+cnr$start
  chrnum_cnr <- as.integer(cnr$chromosome)
  chrnum_cnr[cnr$chromosome=="X"] <- 23
  chrnum_cnr[cnr$chromosome=="Y"] <- 24
  y_cnr <- cnr$log2-chrnum_cnr*8
  
  cns <- read.table(paste("/home/rebecka.bergstrom/cnvkit/1stRunBreast/batchRun/", files_cns[i], sep=""), header = TRUE, as.is=TRUE,sep="\t")
  head(cns)
  chrnum_cns <- as.integer(cns$chromosome)
  chrnum_cns[cns$chromosome=="X"] <- 23
  chrnum_cns[cns$chromosome=="Y"] <- 24
  y_cns <- cns$log2-chrnum_cns*8 
  
  ychrpos <- -(1:24)*8
  title <- paste("log2 copy number ratios ",basename(file_path_sans_ext(files_cnr[i])),sep="")
  filename <- paste("cnvkit/1stRunBreast/Rplots/log2CNratio/",basename(file_path_sans_ext(files_cnr[i])), ".pdf",sep="")
  pdf(filename,width=10,height=10)
  plot(binmiddles,y_cnr, pch=16, col='#00000035',cex=0.6, main=title, yaxt="n", ylab="chromosomes",xlab="position",ylim=c(ychrpos[24]-2,2))
  abline(h=ychrpos,lwd=0.5) #horizontal lines for the chromosomes"zero"-position
  axis(2,at=ychrpos, label=c(1:22,"X","Y")) 
  segments(cns$start,y_cns,cns$end,y_cns,col="#66CD00",lwd=1.8)
  dev.off()
}



