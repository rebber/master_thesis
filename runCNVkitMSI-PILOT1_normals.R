#Run CNVkit for MSI1 normals, use 5 normals as reference
#Created by Rebecka Bergström on 2016-03-15

path_MSI <- "/home/Crisp/clinseq/MSI-PILOT/"
bamFiles_Path <- paste(path_MSI,"*/panel/*B-TD1-CS1-capped.bam",sep="")
ref_samples <- c("ref1B","ref2B","ref3B","ref4B","ref5B")
#exchange "ref1" etc for actual sample IDs
output_dir <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResultsNormals/"
setwd(output_dir)
ref_name <- "reference5MSINormals.cnn"
bamFiles_ref <- c()
for (j in 1:length(ref_samples)) {
  bamFiles_ref[c(2*j-1,2*j)] <- dir(pattern=paste(ref_samples[j],"-TD1-CS1-capped.bam$",sep=""),path_MSI,recursive=TRUE)
}
bamFiles_ref <- bamFiles_ref[-grep("MSI-PILOT",bamFiles_ref)]
bamFiles_ref_Path <- paste(path_MSI,bamFiles_ref," ",sep="",collapse = "")
#systemCommand <- paste("cnvkit.py batch ", bamFiles_Path, " --normal ", bamFiles_ref_Path, "--targets /proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/clinseq_v4.targets.slopped.bed ",
#                       "--split --annotate ~/cnvkit/refFlat_reformatted.txt --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta ",
#                       "--access ~/cnvkit/access-5kb-mappable.hg19_reformatted.bed --output-reference ", ref_name, " ",
#                       "--output-dir ", output_dir, sep="")

systemCommand <- paste("cnvkit.py batch ", bamFiles_Path, " -r /home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResultsNormals/reference5MSINormals.cnn -p 0 --output-dir ", output_dir,sep="")

system(systemCommand)



