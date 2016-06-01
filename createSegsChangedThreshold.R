#Run segmentation in CNVkit with different thresholds for segment breakpoints
#Created by Rebecka Bergström on 2016-03-11

library(tools)

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT2/100kbATbins/full_tumors/"
files_cnr <- dir(pattern="cnr",paste(path_CNVkit,"batchResults",sep=""),recursive=TRUE)
p <- 0.01 #the desired p-value threshold for segmentation
pname <- "p1e-2" #the name of the directory to put the files in
samplenames <- paste(basename(file_path_sans_ext(files_cnr)),"_",pname,".cns",sep="")

for (i in 1:length(files_cnr)) {
  systemCommand <- paste("cnvkit.py segment ", path_CNVkit, "batchResults/",files_cnr[i], " -o ", path_CNVkit, pname, "/cns/",samplenames[i], " -t ",p, sep="")
  system(systemCommand)
}
