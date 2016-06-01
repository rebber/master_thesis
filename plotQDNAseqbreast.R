library(QDNAseq)
library(gplots) 
files <- dir(pattern="T.*qdnaseq.txt", "/home/Crisp/clinseq/BREASTv7/", recursive = TRUE)
savetitle <- c()
for (i in 1:50) {
  dat <- read.table(paste("/home/Crisp/clinseq/BREASTv7/", files[i], sep=""), header = TRUE, as.is=TRUE)
  chr<-dat$chromosome #the chromosome numbers
  chr_ind<-c()
  for (j in 1:22) { #get the indeces where each chromosome starts in chr
    chr_ind[j] <- min(which(chr==j))
  }
  chr_ind <- c(chr_ind, min(which(chr=="X")), min(which(chr=="Y"))) #add indeces for X and Y
  title <- paste0("/home/Crisp/clinseq/BREASTv7/",files[i])
  savetitle[i] <- title
  filename=paste0("breastv7_",i,".pdf")
  pdf(filename)
  plot(dat$segmented,xaxt="n",lwd=0.7,main=title,xlab="Chromosome #",cex.main=0.7)
  abline(v=chr_ind,lwd=0.5) #vertical lines for the chromosomes
  axis(1,at=chr_ind, label=c(1:22,"X","Y"))
  dev.off()
  }

#get the pathways of the files of interest
interesting_files <- c()
j<-0
for (i in c(1,4,6,15,24,27,30,36,42,44)) {
  j<-j+1
  interesting_files[j] <- savetitle [i]
}

# head(dat)
# chr<-dat$chromosome #the chromosome numbers
# max(chr)
# # chr1ind<-which(chr==1)
# # chr2ind<-which(chr==2)
# # chrnum<-c(1:23,"X","Y")
# 
# chr_ind<-c()
# for (i in 1:22) { #get the indeces where each chromosome starts in chr
#   chr_ind[i] <- min(which(chr==i))
# }
# chr_ind <- c(chr_ind, min(which(chr=="X")), min(which(chr=="Y"))) #add indeces for X and Y
# 
#  plot(dat$segmented,xaxt="n",lwd=0.7,main=title,xlab="Chromosome #",cex.main=0.7)
#  abline(v=chr_ind,lwd=0.5) #vertical lines for the chromosomes
#  axis(1,at=chr_ind, label=c(1:22,"X","Y"))
# #axis(2,at=dat$segmented)
