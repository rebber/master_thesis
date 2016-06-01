#Add read count per (CNVkit-)bin to target and antitarget bed-files for panel data
library(tools)

path_bams <- "/home/Crisp/clinseq/MSI-PILOT/"
pattern <- "-TD1-CS1-capped.bam$"
samples <- unique(basename(file_path_sans_ext(dir(pattern=pattern,path_bams,recursive=TRUE)))) #the MSI1-samples
samples <- gsub("-TD1-CS1-capped","",samples)
path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/"

#the bed-files with bins for targets and antitargets
targets <- paste(path_CNVkit,"batchResultsRef5/reference5MSINormals.target-tmp.bed",sep="")
antitargets <- paste(path_CNVkit,"batchResultsRef5/reference5MSINormals.antitarget-tmp.bed",sep="")

for (i in 1:length(samples)) {
  bamFile <- dir(pattern=paste(samples[i],pattern,sep=""),path_bams,recursive=TRUE)
  bamFile <- bamFile[-grep("MSI-PILOT",bamFile)]
  systemCommand_targets <- paste("bedtools intersect -c -a ", targets, " -b ",path_bams, bamFile, " >> ", path_CNVkit, "readCounts/", samples[i],"_readcount.targets.txt", sep="")
  system(systemCommand_targets) #create read count file for targets 
  systemCommand_antitargets <- paste("bedtools intersect -c -a ", antitargets, " -b ",path_bams, bamFile, " >> ", path_CNVkit, "readCounts/", samples[i],"_readcount.antitargets.txt", sep="")
  system(systemCommand_antitargets) #create read count file for antitargets 
}
