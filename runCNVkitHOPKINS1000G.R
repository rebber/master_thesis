#Run CNVkit for HOPKINS samples (panel hopkins_1000G), use reference from 50 normals
#Created by Rebecka Bergström on 2016-05-26

path_main <- "/home/Crisp/clinseq/HOPKINS_PER_KIT/hopkins_1000G/"
output_dir <- "/home/rebecka.bergstrom/cnvkit/HOPKINS/hopkins_1000G/batchResults/"
setwd(output_dir)
ref_name <- paste(output_dir, "reference50normalsPanelHopkins1000G.cnn",sep="")
panel <- "/home/Crisp/clinseq/HOPKINS_PER_KIT/hopkins_1000G_target.bed"

#run CNVkit
systemCommand <- paste("cnvkit.py batch ", path_main, "*.bam --normal ", path_main, "*N.bam --targets ", panel, " --split ",
                       "--annotate ~/cnvkit/refFlat_reformatted.txt --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta ", 
                       "--access ~/cnvkit/access-5kb-mappable.hg19_reformatted.bed --output-reference ", ref_name, " --output-dir ", output_dir, " -p 8", sep="")
system(systemCommand) 



