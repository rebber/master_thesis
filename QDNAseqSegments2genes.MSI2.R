#Create and run the commands for getting log2 CN for genes from QDNAseq segments for wgs in MSI2 samples
#Created by Rebecka Bergströ on 2016-03-10, based on QDNAseqSegments2genes.MSI.R

library(tools)

path_MSI <- "/home/Crisp/clinseq/MSI-PILOT2/"
path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT2/100kbATbins/full_tumors/"

files <- dir(pattern="qdnaseq.segments.txt", path_MSI, recursive = TRUE)
files <- files[-grep("out",files)]
files <- files[-grep("MSI",files)]
samplenames <- dirname(dirname(files))
samplenames <- gsub("-TD3.*-CS1","-wgs",samplenames)

#Command for system call to create bed-file with log2CN-ratio segmented level for each gene, annotated by gene-id (if gene-names are wanted change "-n gene_id" to "-n gene_name")
commands <- paste("python /home/Crisp/rebecka/qdnaseq2bed.py -i ",path_MSI, files, " | sort -k1,1 -k2,2n | 
bedtools map -c 5 -o mean -a /proj/b2010040/private/nobackup/autoseq-genome/genes/ensembl_genes_v75_cleaned_sorted.gtf -b - | 
python /home/Crisp/rebecka/cnvgtf2bed.py -i /dev/stdin -n gene_id > ", path_CNVkit, "segments2genes/", samplenames, ".genes.bed",sep="")

#Call the system
for (i in 1:length(commands)) {
  system(commands[i]) 
}

#Rmember! The segmented values from QDNAseq are not yet log2-transformed! 
