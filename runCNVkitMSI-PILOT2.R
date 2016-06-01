#Run CNVkit for MSI2 samples, use 7 flat tumors as reference, set average bins size for antitargets to 1 Mb
#Created by Rebecka Bergström on 2016-03-08

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

path_MSI <- "/home/Crisp/clinseq/MSI-PILOT2/"
bamFiles_Path <- "/home/Crisp/clinseq/MSI-PILOT2/*/bams/panel/*T-TD3-CB1-panel.bam"
tumRef_samples <- c("ref1T","ref2T","ref3T","ref4T","ref5T","ref6T","ref7T","ref8T","ref9T","ref10T")
#exchange "ref1" etc for actual sample IDs
output_dir <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT2/100kbATbins/full_tumors/batchResults"
setwd(output_dir)
ref_name <- "reference10MSITumors.cnn"
bamFiles_ref <- c()
for (j in 1:length(tumRef_samples)) {
  bamFiles_ref[c(2*j-1,2*j)] <- dir(pattern=paste(tumRef_samples[j],"-TD3-CB1-panel.bam$",sep=""),path_MSI,recursive=TRUE)
}
bamFiles_ref <- bamFiles_ref[-grep("MSI-PILOT2",bamFiles_ref)]
bamFiles_ref_Path <- paste(path_MSI,bamFiles_ref," ",sep="",collapse = "")
systemCommand <- paste("cnvkit.py batch ", bamFiles_Path, " --normal ", bamFiles_ref_Path, "--targets /proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/big_design.targets.slopped.bed ",
                       "--split --annotate ~/cnvkit/refFlat_reformatted.txt --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta ",
                       "--access ~/cnvkit/access-5kb-mappable.hg19_reformatted.bed --output-reference " ref_name, " ",
                       "--output-dir ", output_dir, sep="")

system(systemCommand)


