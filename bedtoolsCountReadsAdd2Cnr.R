#Add read count per (CNVkit-)bin to target and antitarget bed-files for panel data; then add the read count per bin to the .cnr-files from CNVkit
library(tools)

#count reads per bin in bam-files
path_bams <- "/home/Crisp/clinseq/from_rasta/BREASTv6/"
load("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/test50/test50") #loads the testset variable
samples <- testset[,1] #the Bv6-samples
samples <- gsub("_panel_v1","",samples)
path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/ref327norm/"
setwd(paste(path_CNVkit,"batchResults",sep=""))

#the bed-files with bins for targets and antitargets
targets <- paste(path_CNVkit,"batchResults/clinseq_v1_targets.targets.slopped.target.bed",sep="")
antitargets <- paste(path_CNVkit,"batchResults/clinseq_v1_targets.targets.slopped.antitarget.bed",sep="")
#targets <- paste(path_CNVkit,"batchResults/referenceFlatPanelV1.target-tmp.bed",sep="")
#antitargets <- paste(path_CNVkit,"batchresults/clinseq_v1_targets.targets.slopped.antitarget.bed",sep="")

for (i in 1:length(samples)) {
  bamFile <- dir(pattern=paste(testset[i,1],".bam$",sep=""),path=path_bams,recursive=TRUE)
  end <- length(bamFile)
  #if (length(bamFile)==0) {
  #  print(paste("No bam file for sample", i, samples[i]))
  #} else if (length(bamFile)>1) {
  #  print(paste("The number of bamfiles for sample", i, samples[i], "is", length(bamFile)))
  #}
  systemCommand_targets <- paste("bedtools intersect -c -a ", targets, " -b ",path_bams, bamFile[end], " -sorted > ", path_CNVkit, "readCounts/", samples[i],"_readcount.targets.txt", sep="")
  system(systemCommand_targets) #create read count file for targets 
  systemCommand_antitargets <- paste("bedtools intersect -c -a ", antitargets, " -b ",path_bams, bamFile[end], " -sorted > ", path_CNVkit, "readCounts/", samples[i],"_readcount.antitargets.txt", sep="")
  system(systemCommand_antitargets) #create read count file for antitargets 
  print(paste("Sample", i, samples[i], "is counted.", signif(i/length(samples),digits=2)*100, "% done."))
}

############################

#add read count to .cnr-files

#the files for the panel data
files_cnr <- dir(pattern="cnr", paste(path_CNVkit,"batchResults",sep=""), recursive = TRUE)
files_rc_targets <- dir(pattern="T*readcount.targets.txt", paste(path_CNVkit,"readCounts",sep=""),recursive=TRUE)
files_rc_antitargets <- dir(pattern="T*readcount.antitargets.txt", paste(path_CNVkit,"readCounts",sep=""),recursive=TRUE)

for (i in 1:length(files_cnr)) {
  #the panel data with the tumor reference
  dat_cnr <- read.table(paste(path_CNVkit,"batchResults/", files_cnr[i], sep=""), header = TRUE, as.is=TRUE)
  dat_cnr$reads <- rep(NA, length(dat_cnr$chromosome))
  dat_cnr$targetType <- rep(NA, length(dat_cnr$chromosome))
  dat_rc_targets <- read.table(paste(path_CNVkit,"readCounts/", files_rc_targets[i], sep=""), as.is=TRUE)
  colnames(dat_rc_targets) <- c("chromosome","start","end","gene","readCount")
  dat_rc_targets <- unique(dat_rc_targets)
  dat_rc_antitargets <- read.table(paste(path_CNVkit,"readCounts/", files_rc_antitargets[i], sep=""), as.is=TRUE)
  colnames(dat_rc_antitargets) <- c("chromosome","start","end","gene","readCount")
  dat_rc_antitargets <- unique(dat_rc_antitargets)
  dat_rc <- rbind(dat_rc_targets,dat_rc_antitargets)
  #add the raw read count for each bin
  for (j in 1:length(dat_cnr$chromosome)) {
    readsTarget <- dat_rc_targets$readCount[which(dat_cnr$chromosome[j]==dat_rc_targets$chromosome & dat_cnr$start[j]==dat_rc_targets$start)]
    if (length(readsTarget) > 0) {
      dat_cnr$reads[j] <- readsTarget
      dat_cnr$targetType[j] <- "T"
    } else {
      readsAntitarget <- dat_rc_antitargets$readCount[which(dat_cnr$chromosome[j]==dat_rc_antitargets$chromosome & dat_cnr$start[j]==dat_rc_antitargets$start)]
      if (length(readsAntitarget) > 0) {
        dat_cnr$reads[j] <- readsAntitarget
        dat_cnr$targetType[j] <- "AT"
      } else {
        print(paste("Bin", j, "in sample", files_cnr[i], "had no read count or target type assigned"))
      }
    }
  }
  
  
  #Add smoothed bin values
  dat_cnr$smooth <- rep(NA,length(dat_cnr$chromosome))
  targets <- which(dat_cnr$targetType=="T")
  antitargets <- which(dat_cnr$targetType=="AT")
  for (j in 1:length(dat_cnr$smooth)) { 
    if (dat_cnr$targetType[j]=="T") {
      pos <- which(targets==j)
      if (pos>5 && (length(targets) - pos) > 5) {
        start <- pos - 5
        end <- pos + 5
      } else if (length(targets) - pos <= 5) {
        start <- length(targets) - 11
        end <- length(targets)
      } else if (pos<=5) {
        start <- 1 
        end <- 11
      }
      ind <- targets[start:end]
    } else {
      pos <- which(antitargets==j)
      if (pos>5 && (length(antitargets) - pos) > 5) {
        start <- pos - 5
        end <- pos + 5
      } else if (length(antitargets) - pos <= 5) {
        start <- length(antitargets) - 11
        end <- length(antitargets)
      } else if (pos<=5) {
        start <- 1 
        end <- 11
      }
      ind <- antitargets[start:end]
    }
    dat_cnr$smooth[j] <- median(dat_cnr$log2[ind])
  }
  
  #write a new file which is the .cnr file with raw read count target type and smoothed bin log2CN added
  write.table(dat_cnr,paste(path_CNVkit,"cnrRcComb/addedSmooth/",file_path_sans_ext(files_cnr[i]),".txt",sep=""),quote=FALSE,row.names=FALSE)
  
  print(paste("Counts and smoothed for sample", i, files_cnr[i], "are added.", signif(i/length(files_cnr),digits=2)*100, "% done."))
}


#test_target_bed <- read.table("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/batchResults/clinseq_v1_targets.targets.slopped.target.bed",sep="\t")
#test_target_bed <- read.table("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/batchResults/clinseq_v1_targets.targets.slopped.target.bed",sep=",")
