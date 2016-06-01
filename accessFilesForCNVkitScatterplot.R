#Run cnvkit.py scatter for multiple files

setwd("~/cnvkit/MSI-PILOT/batchResults/")
library(tools)

files_cnr <- dir(pattern=".cnr", "~/cnvkit/MSI-PILOT/batchResults/", recursive = TRUE)
files_cnr <- files_cnr[-grep("broken_pipe",files_cnr)]
files_cns <- dir(pattern=".cns", "~/cnvkit/MSI-PILOT/batchResults/", recursive = TRUE)
files_cns <- files_cns[-grep("broken_pipe",files_cns)]
samplenames <- file_path_sans_ext(files_cnr)
scatterfilenames <- paste(samplenames,"-scatter.ymin.pdf",sep="")

commands <- paste("cnvkit.py scatter", files_cnr, "-s", files_cns, "--y-min -3 -o", scatterfilenames)
#system(command) do not recognize command cnvkit.py, run each row of commands from the command line instead
