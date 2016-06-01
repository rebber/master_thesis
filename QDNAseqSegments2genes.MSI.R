#Create and run the commands for getting log2 CN for genes from QDNAseq segments for wgs in MSI samples

library(tools)

files <- dir(pattern=".*T-TD1-WGS.qdnaseq.txt", "/home/Crisp/clinseq/MSI-PILOT/", recursive = TRUE)
files <- files[-grep("out",files)]
files <- files[-grep("MSI",files)]
samplenames <- basename(file_path_sans_ext(files))


commands <- paste("python /home/Crisp/rebecka/qdnaseq2bed.py -i /home/Crisp/clinseq/MSI-PILOT/", files, " | sort -k1,1 -k2,2n | 
bedtools map -c 5 -o mean -a /proj/b2010040/private/nobackup/autoseq-genome/genes/ensembl_genes_v75_cleaned_sorted.gtf -b - | 
python /home/Crisp/rebecka/cnvgtf2bed.py -i /dev/stdin -n gene_name > /home/rebecka.bergstrom/cnvkit/MSI-PILOT/QDNAseqWGS/segments2genes/", samplenames, ".genes.bed",sep="")

for (i in 1:length(commands)) {
  system(commands[i]) 
}

write.table(commands,"~/cnvkit/scripts/QDNAseqSegments2genes.MSI.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
#I add description and empty rows between each commando with nano
