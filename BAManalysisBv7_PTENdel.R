#Analyze the BAM-files (number of reads etc) for the tumor and normal samples which showed homozygos deletion of PTEN for tumors.

#source("http://bioconductor.org/biocLite.R")
#biocLite("Rsamtools")
library(Rsamtools)
library(tools)

#setwd("/home/rebecka.bergstrom/BREASTv7_PTEN_del/panelBAMfileAnalysis/")
setwd("/home/Crisp/rebecka/breastv7_bams/")

samples <- c("PTENdel1","PTENdel2","PTENdel3","PTENdel4")
#exchange "PTENdel1" etc for correct sample IDs
bamStats <- as.data.frame(cbind(samples,matrix(rep(NA,length(samples)*6),ncol=6)),stringsAsFactors = FALSE)
colnames(bamStats) <- c("sample","TtotReads","TmappedChr","TmappedRel","NtotReads","NmappedChr","NmappedRel")

#BAM-statistics for tumor BAMs
for (i in 1:length(samples)) {
  systemCommand <- paste("samtools idxstats ",samples[i],"T_panel_v1.bam",sep="")
  idxstatsOutput <- data.frame(do.call(rbind,strsplit(system(systemCommand, intern=TRUE),"\t")),stringsAsFactors = FALSE )
  colnames(idxstatsOutput) <- c("sequence","length","mapped","unmapped")
  idxstatsOutput$length <- as.double(idxstatsOutput$length)
  idxstatsOutput$mapped <- as.double(idxstatsOutput$mapped)
  idxstatsOutput$unmapped <- as.double(idxstatsOutput$unmapped)
  write.table(idxstatsOutput, paste("/home/rebecka.bergstrom/BREASTv7_PTEN_del/panelBAMfileAnalysis/",samples[i],"T_panel_v1.idxstats.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)
  
  mapped_tot <- sum(idxstatsOutput$mapped)
  mapped_chr <- sum(idxstatsOutput$mapped[1:24])
  unmapped_tot <- sum(idxstatsOutput$unmapped)
  reads_tot <- mapped_tot+unmapped_tot
  mapped_rel <- mapped_chr/reads_tot
  
  bamStats$TtotReads[i] <- reads_tot
  bamStats$TmappedChr[i] <- mapped_chr
  bamStats$TmappedRel[i] <- mapped_rel
}

#BAM-statistics for normal BAMs
for (i in 1:length(samples)) {
  systemCommand <- paste("samtools idxstats ",samples[i],"N_panel_v1.bam",sep="")
  idxstatsOutput <- data.frame(do.call(rbind,strsplit(system(systemCommand, intern=TRUE),"\t")),stringsAsFactors = FALSE )
  colnames(idxstatsOutput) <- c("sequence","length","mapped","unmapped")
  idxstatsOutput$length <- as.double(idxstatsOutput$length)
  idxstatsOutput$mapped <- as.double(idxstatsOutput$mapped)
  idxstatsOutput$unmapped <- as.double(idxstatsOutput$unmapped)
  write.table(idxstatsOutput, paste("/home/rebecka.bergstrom/BREASTv7_PTEN_del/panelBAMfileAnalysis/",samples[i],"N_panel_v1.idxstats.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)
  
  mapped_tot <- sum(idxstatsOutput$mapped)
  mapped_chr <- sum(idxstatsOutput$mapped[1:24])
  unmapped_tot <- sum(idxstatsOutput$unmapped)
  reads_tot <- mapped_tot+unmapped_tot
  mapped_rel <- mapped_chr/reads_tot
  
  bamStats$NtotReads[i] <- reads_tot
  bamStats$NmappedChr[i] <- mapped_chr
  bamStats$NmappedRel[i] <- mapped_rel
}

#Change variable class from character to double for the statistics in bamStats so it can be used
bamStats$TtotReads <- as.double(bamStats$TtotReads)
bamStats$TmappedChr <- as.double(bamStats$TmappedChr)
bamStats$TmappedRel <- as.double(bamStats$TmappedRel)
bamStats$NtotReads <- as.double(bamStats$NtotReads)
bamStats$NmappedChr <- as.double(bamStats$NmappedChr)
bamStats$NmappedRel <- as.double(bamStats$NmappedRel)



#Dilute the tumor-reads with reads from matched normal
#tumProp <- 0.75 #the proportion of tumor reads we want in the resulting BAM file after dilution
tumPerc <- 75 #the percentage of tumor reads we want in the resulting BAM file after dilution

for (i in 2:length(samples)) {
  fileT <- paste("/home/Crisp/rebecka/breastv7_bams/",samples[i],"T_panel_v1.bam",sep="")
  fileN <- paste("/home/Crisp/rebecka/breastv7_bams/",samples[i],"N_panel_v1.bam",sep="")
  fileNpart <- paste("/home/rebecka.bergstrom/BREASTv7_PTEN_del/panelBAMfileAnalysis/dilutions/tumorPerc_",tumPerc,"/",samples[i],"N_partOfReads.bam",sep="")
  fileDil <- paste("/home/rebecka.bergstrom/BREASTv7_PTEN_del/panelBAMfileAnalysis/dilutions/tumorPerc_",tumPerc,"/",samples[i],"_diluted_",tumPerc,".bam",sep="")
  
  numNormReads <- (100-tumPerc)*bamStats$TtotReads[i]/tumPerc #the number of reads to use from the normal BAM
  normProp <- numNormReads/bamStats$NtotReads[i] #the proportion of reads to extract from normal BAM
  systemCommand <- paste("samtools view -s ",normProp," -b ", fileN," > ",fileNpart, sep="") #paste command to the system
  sysOutNormProp <- system(systemCommand,intern=TRUE) #call the command line with the samtools view command to create a new BAM file with the reads from the normal which should be used for the dilution
  mergeBam(c(fileT,fileNpart),fileDil, indexDestination=TRUE)
  
  #check the created, diluted bam-file
  systemCommand <- paste("samtools idxstats ",fileDil,sep="")
  idxstatsOutputDiluted <- data.frame(do.call(rbind,strsplit(system(systemCommand, intern=TRUE),"\t")),stringsAsFactors = FALSE )
  colnames(idxstatsOutputDiluted) <- c("sequence","length","mapped","unmapped")
  idxstatsOutputDiluted$length <- as.double(idxstatsOutputDiluted$length)
  idxstatsOutputDiluted$mapped <- as.double(idxstatsOutputDiluted$mapped)
  idxstatsOutputDiluted$unmapped <- as.double(idxstatsOutputDiluted$unmapped)
  
  sum(idxstatsOutputDiluted$mapped)
  sum(idxstatsOutputDiluted$unmapped)
  sum(idxstatsOutputDiluted$mapped)+sum(idxstatsOutputDiluted$unmapped)
  
  #Is the #tumor reads 75 % of all reads in diluted file? Yes, the portion of tumor reads represents 75.00274 %. 
  bamStats$TtotReads[i]/(sum(idxstatsOutputDiluted$mapped)+sum(idxstatsOutputDiluted$unmapped)) 
}
