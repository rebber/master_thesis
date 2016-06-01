#An attempt to run CNVkit from an R interface, to avoid switching between R and command line, and enable scripts for the whole workflow in dilution with bam-analysis, CNVkit, plotting, t-testing etc.

#Add miniconda2/bin to the environment paths to enable CNVkit call via system()
presentPath <- Sys.getenv("PATH")
cnvkitAddedPath <- paste("/home/rebecka.bergstrom/miniconda2/bin:", presentPath,sep="",collapse="")
Sys.setenv(PATH=cnvkitAddedPath)
#Add the R library in miniconda which CNVkit uses, to the libraries which this R environment searches
presentRLibUser <- Sys.getenv("R_LIBS_USER")
presentRLibSites <- Sys.getenv("R_LIBS_SITE")
minicondaRAddedPath <- paste("/home/rebecka.bergstrom/miniconda2/lib/R/library:",presentRLibUser,sep="")
Sys.setenv(R_LIBS_USER=minicondaRAddedPath)
Sys.setenv(R_LIBS_USER="/home/rebecka.bergstrom/miniconda2/lib/R/library")
Sys.setenv(R_HOME="/home/rebecka.bergstrom/miniconda2/lib/R")


#tumProp <- 0.75 #the proportion of tumor reads in diluted bam files
tumPerc <- 75 #the percentage of tumor reads in diluted bam files
bamFiles <- dir(pattern=paste("diluted_",tumPerc,sep=""),paste("/home/rebecka.bergstrom/BREASTv7_PTEN_del/panelBAMfileAnalysis/dilutions/tumorPerc_",tumPerc,sep=""),recursive=TRUE)
bamFiles <- bamFiles[-grep("bai",bamFiles)]
bamFilesPath <- paste("/home/rebecka.bergstrom/BREASTv7_PTEN_del/panelBAMfileAnalysis/dilutions/tumorPerc_", tumPerc, "/",bamFiles," ",sep="",collapse="")
systemCommand <- paste("cnvkit.py batch ", bamFilesPath, "-r /home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/reference14breast.cnn -d /home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/dilutions/tumorPerc_",tumPerc, "/batchResults/",sep="")
sysOutCnvkit <- system(systemCommand, intern = TRUE)


