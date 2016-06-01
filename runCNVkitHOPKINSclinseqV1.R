#Run CNVkit for HOPKINS samples (panel clinseq_v1), use reference from 50 normals
#Created by Rebecka Bergström on 2016-05-26

path_main <- "/home/Crisp/clinseq/HOPKINS_PER_KIT/clinseq_v1/"
output_dir <- "/home/rebecka.bergstrom/cnvkit/HOPKINS/clinseq_v1/batchResults/"
setwd(output_dir)
ref_name <- paste(output_dir, "reference50normalsPanelV1.cnn",sep="")
panel <- "/home/Crisp/clinseq/HOPKINS_PER_KIT/clinseq_v1_targets.bed"

#run CNVkit
systemCommand <- paste("cnvkit.py batch ", path_main, "*.bam --normal ", path_main, "*N.bam --targets ", panel, " --split ",
                       "--annotate ~/cnvkit/refFlat_reformatted.txt --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta ", 
                       "--access ~/cnvkit/access-5kb-mappable.hg19_reformatted.bed --output-reference ", ref_name, " --output-dir ", output_dir, " -p 8", sep="")
system(systemCommand) 



