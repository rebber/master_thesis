#Run CNVkit for BREASTv6 tumors, use flat reference
#Created by Rebecka Bergström on 2016-03-22


path_main <- "/home/Crisp/clinseq/from_rasta/BREASTv6/"
bamFiles_Path <- paste(path_main,"*/datareports/*/*T_panel_v1.bam",sep="") # tumor files
output_dir <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/batchResults/"
setwd(output_dir)
ref_name <- paste(output_dir, "referenceFlatPanelV1.cnn",sep="")
panel <- "/proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/clinseq_v1_targets.targets.slopped.bed"
targets <- paste(output_dir,"clinseq_v1_targets.targets.slopped.target.bed",sep="")
antitargets <- paste(output_dir,"clinseq_v1_targets.targets.slopped.antitarget.bed",sep="")

#build the reference
systemCommand <- paste("cnvkit.py target ", panel, " --annotate ~/cnvkit/refFlat_reformatted.txt --split -o ", targets, sep="")
system(systemCommand)
systemCommand <- paste("cnvkit.py antitarget ", targets, " --access ~/cnvkit/access-5kb-mappable.hg19_reformatted.bed -o ", antitargets,sep="")
system(systemCommand)
systemCommand <- paste("cnvkit.py reference -o ", ref_name, " -t ", targets, " -a ", antitargets, " --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta",sep="")
system(systemCommand)

#run the tumors 
systemCommand <- paste("cnvkit.py batch ", bamFiles_Path, " -r ", ref_name, " -p 5 --output-dir ", output_dir,sep="")
system(systemCommand)



