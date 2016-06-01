#Create and run the commands for getting log2 CN for genes from CNVkitsegments for panel in MSI samples

library(tools)

files_cns <- dir(pattern="cns", "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResults/", recursive = TRUE)
files_cns <- files_cns[-grep("broken_pipe",files_cns)]
samplenames_cns <- basename(file_path_sans_ext(files_cns))

for (i in 1:length(files_cns)) {
  data_cns <- read.table(paste("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResults/",files_cns[i],sep=""), as.is=TRUE,header=TRUE)
  colnames(data_cns) <- c("chromosome","start","end","gene","segmented","probes","weight")
  write.table(data_cns,paste("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/segments2genes/",samplenames_cns[i],".txt",sep=""),row.names=FALSE,quote=FALSE,sep="\t")
}

files <-dir(pattern="txt","/home/rebecka.bergstrom/cnvkit/MSI-PILOT/segments2genes/",recursive=TRUE)
samplenames <- basename(file_path_sans_ext(files))


commands <- paste("python /home/Crisp/rebecka/qdnaseq2bed.py -i /home/rebecka.bergstrom/cnvkit/MSI-PILOT/segments2genes/", files, " | sort -k1,1 -k2,2n | 
bedtools map -c 5 -o mean -a /proj/b2010040/private/nobackup/autoseq-genome/genes/ensembl_genes_v75_cleaned_sorted.gtf -b - | 
python /home/Crisp/rebecka/cnvgtf2bed.py -i /dev/stdin -n gene_name > /home/rebecka.bergstrom/cnvkit/MSI-PILOT/segments2genes/", samplenames, ".genes.bed",sep="")

for (i in 1:length(commands)) {
  system(commands[i]) 
}

write.table(commands,"~/cnvkit/scripts/CNVkitSegments2genes.MSI.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
#I add description and empty rows between each commando with nano