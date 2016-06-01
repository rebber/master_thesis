# Create files that contain the CNVkit .cnr data, the raw read count for each bin and an annotation of which bin type it is
#Part 4 of 4, file 211-281
#Created by Rebecka Bergström on 2016-03-23

library(tools)
path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/"

#the files for the panel data
files_cnr <- dir(pattern="cnr", paste(path_CNVkit,"batchResults",sep=""), recursive = TRUE)
files_rc_targets <- dir(pattern="T_readcount.targets.txt", paste(path_CNVkit,"readCounts",sep=""),recursive=TRUE)
files_rc_antitargets <- dir(pattern="T_readcount.antitargets.txt", paste(path_CNVkit,"readCounts",sep=""),recursive=TRUE)

for (i in 211:281) {
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
  
  #write a new file which is the .cnr file with raw read count and target type added
  write.table(dat_cnr,paste(path_CNVkit,"cnrRcComb/",file_path_sans_ext(files_cnr[i]),".txt",sep=""),quote=FALSE,row.names=FALSE)
  
  print(paste("Counts for sample", i, files_cnr[i], "are added.", 100*signif(i/length(files_cnr),digits=2), "% done."))
}