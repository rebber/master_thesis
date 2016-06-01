#Run CNVkit for MSI1 tumors and normals, use flat reference
#Created by Rebecka Bergström on 2016-03-16

path_MSI <- "/home/Crisp/clinseq/MSI-PILOT/"
bamFiles_Path <- paste(path_MSI,"*/panel/*-TD1-CS1-capped.bam",sep="") #both normal and tumor files
output_dir <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResultsRefFlat/"
setwd(output_dir)
ref_name <- "referenceFlat.cnn"
targets <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResultsRef5/reference5MSINormals.target-tmp.bed"
antitargets <- "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/batchResultsRef5/reference5MSINormals.antitarget-tmp.bed"

#build the reference
systemCommand <- paste("cnvkit.py reference -o ", ref_name, " -t ", targets, " -a ", antitargets, " --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta ",sep="")
system(systemCommand)

#run the tumors and normals
systemCommand <- paste("cnvkit.py batch ", bamFiles_Path, " -r ",output_dir, ref_name, " -p 0 --output-dir ", output_dir,sep="")
system(systemCommand)



