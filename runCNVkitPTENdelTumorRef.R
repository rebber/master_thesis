#Run CNVkit for Bv7 samples with homoz PTEN deletion, use 11 flat tumors as reference

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
bamFiles_PTENdel <- c()
for (i in 1:length(samples)) {
  bamFiles_PTENdel <- append(bamFiles_PTENdel,dir(pattern=paste(samples[i],"T_panel_v1.bam",sep=""),"/home/Crisp/rebecka/breastv7_bams/",recursive=TRUE))
}
bamFiles_ref <- dir(pattern="T.*panel_v1.bam","/home/Crisp/rebecka/breastv7_bams_with_no_cna/",recursive=TRUE)
bamFiles_PTENdel_Path <- paste("/home/Crisp/rebecka/breastv7_bams/",bamFiles_PTENdel," ",sep="",collapse="")
bamFiles_ref_Path <- paste("/home/Crisp/rebecka/breastv7_bams_with_no_cna/",bamFiles_ref," ",sep="",collapse = "")
systemCommand <- paste("cnvkit.py batch ", bamFiles_PTENdel_Path, "--normal ", bamFiles_ref_Path, "--targets /proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/clinseq_v1_targets.targets.slopped.bed ",
                       "--split --annotate ~/cnvkit/refFlat_reformatted.txt --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta ",
                       "--access ~/cnvkit/access-5kb-mappable.hg19_reformatted.bed --output-reference reference11breastTumors.cnn ",
                       "--output-dir /home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/tumRef/full_tumors/batchResults/",sep="")
sysOutCnvkit <- system(systemCommand, intern = TRUE)


