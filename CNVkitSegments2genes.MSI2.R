#Create and run the commands for getting log2 CN for genes from CNVkitsegments for panel in MSI2 samples
#Created by Rebecka Bergströ on 2016-03-10, based on CNVkitSegments2genes.MSI.R

library(tools)

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT2/100kbATbins/full_tumors/"
files_cns <- dir(pattern="cns", paste(path_CNVkit,"batchResults/",sep=""), recursive = TRUE)
samplenames_cns <- basename(file_path_sans_ext(files_cns))

#Create txt-files with appropriate headers (qdnaseq2bed.py takes txt-file with one column "segmented") from the cns-files from CNVkit
#Remember!  These segmented values are already log2-transformed!
for (i in 1:length(files_cns)) {
  data_cns <- read.table(paste(path_CNVkit,"batchResults/",files_cns[i],sep=""), as.is=TRUE,header=TRUE)
  colnames(data_cns) <- c("chromosome","start","end","gene","segmented","probes","weight")
  write.table(data_cns,paste(path_CNVkit,"segments2genes/",samplenames_cns[i],".txt",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
}

files <-dir(pattern="txt",paste(path_CNVkit,"segments2genes/",sep=""),recursive=TRUE)
samplenames <- basename(file_path_sans_ext(files))

#Command for system call to create bed-file with log2CN-ratio segmented level for each gene, annotated by gene-id (if gene-names are wanted change "-n gene_id" to "-n gene_name")
commands <- paste("python /home/Crisp/rebecka/qdnaseq2bed.py -i ", path_CNVkit, "segments2genes/", files, " | sort -k1,1 -k2,2n | 
bedtools map -c 5 -o mean -a /proj/b2010040/private/nobackup/autoseq-genome/genes/ensembl_genes_v75_cleaned_sorted.gtf -b - | 
python /home/Crisp/rebecka/cnvgtf2bed.py -i /dev/stdin -n gene_id > ", path_CNVkit, "segments2genes/", samplenames, ".genes.bed",sep="")

#Call the system
for (i in 1:length(commands)) {
  system(commands[i]) 
}


