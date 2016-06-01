
setwd("/home/Crisp/rebecka/")
library(tools)
files <- dir(pattern="T.*qdnaseq.txt", "/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", recursive = TRUE)
for (i in 1:length(files)) {
  dat <- read.table(paste("/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", files[i], sep=""), header = TRUE, as.is=TRUE)
  log2CN <- log2(dat$copynumber)
  position <- dat$start
  chr<-dat$chromosome #the chromosome numbers
  #chr_ind<-c()
  #for (j in 1:22) { #get the indeces where each chromosome starts in chr
  #  chr_ind[j] <- min(which(chr==j))
  #}
  #chr_ind <- c(chr_ind, min(which(chr=="X")), min(which(chr=="Y")), max(which(chr=="Y"))) #add indeces for beginning of X and Y and end of Y
  #title <- paste0("/home/Crisp/clinseq/BREASTv7/",files[i])
  ychr <- -as.integer(chr)*8
  ychrpos <- -(1:22)*8
  title <- paste("log2 copy number ",basename(file_path_sans_ext(files[i])),sep="")
  filename <- paste(basename(file_path_sans_ext(files[i])), ".pdf",sep="")
  pdf(filename,width=10,height=10)
  #plot(log2CN,xaxt="n",lwd=0.7,main=title,xlab="Chromosome #" ) #Om titeln ska vara hela sökvägen sätt cex.main=0.7 så att den får plats
  if (grep("wgs",filename)==1){
    plot(position,log2CN+ychr, pch=16, col='#00000003',cex=0.4, main=title, yaxt="n", ylab="chromosomes")
    abline(h=ychrpos,lwd=0.5) #horizontal lines for the chromosomes"zero"-position
    axis(2,at=ychrpos, label=1:22) 
  }
  else {
    plot(position,log2CN+ychr, pch=16, col='#00000020',cex=0.5, main=title, yaxt="n", ylab="chromosomes")
    abline(h=ychrpos,lwd=0.5) #horizontal lines for the chromosomes"zero"-position
    axis(2,at=ychrpos, label=1:22) 
  }
  dev.off()
}


