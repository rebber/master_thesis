#Find the samples from experiment BREASTv7 that have an amplification in ERBB2

#Find which bins (from QDNAseq data, 15Kb bins) overlaps with ERBB2
genes_bed <- read.table("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/QDNAseqWGS/segments2genes/example_file.qdnaseq.genes.bed",as.is=TRUE) #Use annotation from another sample (the annotation is the same, ignore CN data)
ERBB2_ind <- grep("ERBB2",genes_bed$V4) #Finds both ERBB2 and ERBB2IP
ERBB2_ind <- ERBB2_ind[-which(ERBB2_ind==grep("ERBB2IP",genes_bed$V4))] #Removes the index for ERBB2IP
ERBB2_bed <- genes_bed[ERBB2_ind,] #The bed-data for ERBB2
colnames(ERBB2_bed) <- c("chromosome", "start", "end","gene","CN","strand") #the column CN is to ignore since this annotation is gotten from another sample (MSI-PILOT), but the gene annotation and localization is the same

data_1stSample <- read.table("/home/Crisp/clinseq/BREASTv7/a1/a1Textr2/wgs/a1Textr2_wgs.qdnaseq.txt",as.is=TRUE,header=TRUE) #Use the first sample in this cohort to find the correct bin division
#exchange "a1" for correct sample IDs
data_1stSample_chr17 <- data_1stSample[which(data_1stSample$chromosome==ERBB2_bed$chromosome),] #ERBB2 is on chr 17
ERBB2_startcoordChr <- floor(ERBB2_bed$start/15000)*15000+1 #the start pos on the chrom of the 1st bin overlapping with ERBB2
ERBB2_endcoordChr <- ceiling(ERBB2_bed$end/15000)*15000 #the end pos on the chrom of the last bin overlapping with ERBB2
ERBB2_startind <- which(data_1stSample$chromosome==ERBB2_bed$chromosome & data_1stSample$start==ERBB2_startcoordChr) #the overall index of the 1st bin overlapping w/ ERBB2
ERBB2_endind <- which(data_1stSample$chromosome==ERBB2_bed$chromosome & data_1stSample$end==ERBB2_endcoordChr) #the overall index of the last bin overlapping w/ ERBB2
numBins_ERBB2 <- ERBB2_endind-ERBB2_startind+1 #Four bins overlap in this case

#Search the samples in BREASTv7
files <- dir(pattern="T.*wgs.qdnaseq.txt","/home/Crisp/clinseq/BREASTv7/",recursive=TRUE)
ERBB2_CNs <- as.data.frame(cbind(files,matrix(rep(NA,length(files)*numBins_ERBB2),ncol=numBins_ERBB2)),stringsAsFactors=FALSE)

for (i in 1:length(files)) {
  dat <- read.table(paste("/home/Crisp/clinseq/BREASTv7/",files[i],sep=""),as.is=TRUE,header=TRUE)
  ERBB2_CNs[i,2:(numBins_ERBB2+1)] <- dat$copynumber[ERBB2_startind:ERBB2_endind]
}

ERBB2_medians <- cbind(files,rep(NA,length(files)))
for (i in 1:length(ERBB2_CNs$files)) {
  ERBB2_medians[i,2] <- median(as.numeric(ERBB2_CNs[i,2:(numBins_ERBB2+1)]))
}

quantile(as.numeric(as.matrix(ERBB2_CNs[,2:(numBins_ERBB2+1)])),probs=seq(0,1,0.125))
quantile(as.numeric(ERBB2_medians[,2]))
sum(as.numeric(ERBB2_medians[,2])>10)
sum(log2(as.numeric(ERBB2_medians[,2]))>3)

ERBB2_amplification <- ERBB2_medians[(which(log2(as.numeric(ERBB2_medians[,2]))>3)),1]
write.table(ERBB2_amplification, "/home/rebecka.bergstrom/BREASTv7_PTEN_del/ERBB2amp/BREASTv7_ERBB2_log2CN_over_3.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)



