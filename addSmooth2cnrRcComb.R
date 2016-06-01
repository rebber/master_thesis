#Add "smoothed" bin value, median based on 10 closest bins of the same type for each bin
#Created by Rebecka Bergström on 2016-04-07

library(tools)

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/"

#find the files
files_cnr <- dir(pattern="T*panel_v1.txt", paste(path_CNVkit,"cnrRcComb",sep=""), recursive = TRUE)
files_cnr <- files_cnr[-grep("T_panel",files_cnr)]

for (i in 1:length(files_cnr_tumors_refFlat)) { 
  ptm <- proc.time()
  #with flat reference:
  #the data for the tumors_refFlat
  dat_cnr <- read.table(paste(path_CNVkit,"cnrRcComb/", files_cnr[i], sep=""), header = TRUE, as.is=TRUE)
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
  write.table(dat_cnr,paste(path_CNVkit,"cnrRcComb/addedSmooth/",file_path_sans_ext(files_cnr[i]),".txt",sep=""),quote=FALSE,row.names=FALSE)
  print(i)
  print(proc.time()-ptm)
}  