#Create new .cnr files where log2-values are exchanged to pca-adjusted log2-values.
#Created by Rebecka Bergström on 2016-04-13

library(tools)

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/"
load(paste(path_CNVkit,"test50/test50",sep="")) #load the variable testset, containing names and indeces of test samples
files_cnr_all <- dir(pattern=".cnr",paste(path_CNVkit,"allCnr/",sep=""),recursive=TRUE)
files_cnr <- files_cnr_all[as.integer(testset[,2])]

adj_log2 <- as.matrix(read.table(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/adj_log2_testset.txt",sep=""),as.is=TRUE))

for (i in 1:length(files_cnr)) {
  dat_cnr <- read.table(paste(path_CNVkit,"allCnr/",files_cnr[i],sep=""),header=TRUE,as.is = TRUE)
  Ybins <- which(dat_cnr$chromosome=="Y") 
  dat_cnr$log2 <- c(adj_log2[i,],dat_cnr$log2[Ybins]) #change to adjusted log2-values in all bins except on Y chrom where such are not avilable
  write.table(dat_cnr,file=paste(path_CNVkit,"segmAdjBins/segmResults/",files_cnr[i],sep=""),quote=FALSE,sep="\t",row.names=FALSE)
}