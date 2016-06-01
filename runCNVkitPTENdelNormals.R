#Run CNVkit for normals for the samples with homozygous PTEN deletion

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


samples <- c("PTENdel1","PTENdel2","PTENdel3","PTENdel4")
#exchange "PTENdel1" etc towards correct sample IDs

bamFiles <- c()
for (i in 1:length(samples)) {
  bamFiles <- append(bamFiles,dir(pattern=paste(samples[i],"N_panel_v1.bam",sep=""),"/home/Crisp/rebecka/breastv7_bams/",recursive=TRUE))
}
#bamFiles <- bamFiles[-grep("bai",bamFiles)]
bamFilesPath <- paste("/home/Crisp/rebecka/breastv7_bams/",bamFiles," ",sep="",collapse="")
systemCommand <- paste("cnvkit.py batch ", bamFilesPath, "-r /home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/reference14breast.cnn -d /home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/batchResults/",sep="")
sysOutCnvkit <- system(systemCommand, intern = TRUE)


