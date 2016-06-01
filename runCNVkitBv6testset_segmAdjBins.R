#Use CNVkit to redo segmentation after bin-values have been adjusted with pca. 
#Created by Rebecka Bergström on 2016-04-13

library(tools)

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/segmAdjBins/segmResults/"
files_cnr <- dir(pattern=".cnr",path=path_CNVkit,recursive=TRUE) #the .cnr files with adjusted log2 bin values
files_cns <- paste(basename(file_path_sans_ext(files_cnr)), ".cns", sep="")
setwd(path_CNVkit)

for (i in 1:length(files_cnr)) {
  systemCommand <- paste("cnvkit.py segment", files_cnr[i], "-o", files_cns[i])
  system(systemCommand)
}