#Find the samples from experiment BREASTv7 that have a deletion in PTEN

#Find which bins (from QDNAseq data, 15Kb bins) overlaps with PTEN
genes_bed <- read.table("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/QDNAseqWGS/segments2genes/example_file.qdnaseq.genes.bed",as.is=TRUE)
PTEN_ind <- grep("PTEN",genes_bed$V4) #Finds both PTEN and PTENP1
PTEN_ind <- PTEN_ind[-which(PTEN_ind==grep("PTENP1",genes_bed$V4))] 
PTEN_bed <- genes_bed[PTEN_ind,] #The bed-data for PTEN
colnames(PTEN_bed) <- c("chromosome", "start", "end","gene","CN","strand")

data_1stSample <- read.table("/home/Crisp/clinseq/BREASTv7/a1/a1Textr2/wgs/a1Textr2_wgs.qdnaseq.txt",as.is=TRUE,header=TRUE) #exchange a1 for correct sample IDs
data_1stSample_chr10 <- data_1stSample[which(data_1stSample$chromosome==PTEN_bed$chromosome),]
PTEN_startcoordChr <- floor(PTEN_bed$start/15000)*15000+1 #the start pos on the chrom of the 1st bin overlapping with PTEN
PTEN_endcoordChr <- ceiling(PTEN_bed$end/15000)*15000 #the end pos on the chrom of the last bin overlapping with PTEN
PTEN_startind <- which(data_1stSample$chromosome==PTEN_bed$chromosome & data_1stSample$start==PTEN_startcoordChr) #the overall index of the 1st bin overlapping w/ PTEN
PTEN_endind <- which(data_1stSample$chromosome==PTEN_bed$chromosome & data_1stSample$end==PTEN_endcoordChr) #the overall index of the last bin overlapping w/ PTEN
numBins_PTEN <- PTEN_endind-PTEN_startind+1 #Nine bins overlap in this case

#Search the samples in BREASTv7
files <- dir(pattern="T.*wgs.qdnaseq.txt","/home/Crisp/clinseq/BREASTv7/",recursive=TRUE)
PTEN_CNs <- as.data.frame(cbind(files,matrix(rep(NA,length(files)*numBins_PTEN),ncol=numBins_PTEN)),stringsAsFactors=FALSE)

for (i in 1:length(files)) {
  dat <- read.table(paste("/home/Crisp/clinseq/BREASTv7/",files[i],sep=""),as.is=TRUE,header=TRUE)
  PTEN_CNs[i,2:(numBins_PTEN+1)] <- dat$copynumber[PTEN_startind:PTEN_endind]
}

PTEN_medians <- cbind(files,rep(NA,length(files)))
for (i in 1:length(PTEN_CNs$files)) {
  PTEN_medians[i,2] <- median(as.numeric(PTEN_CNs[i,2:(numBins_PTEN+1)]))
}

quantile(as.numeric(as.matrix(PTEN_CNs[,2:10])))
quantile(as.numeric(PTEN_medians[,2]))
sum(as.numeric(PTEN_medians[,2])<0.75)

PTEN_deletion <- PTEN_medians[(which(PTEN_medians[,2]<0.75)),1]
write.table(PTEN_deletion, "/home/rebecka.bergstrom/BREASTv7_PTEN_del/BREASTv7_PTEN_CN_under_075.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

log2(0.75)
log2(0.8)

