#Run CNVkit for BREASTv6 tumors with deviant names, use flat reference
#Created by Rebecka Bergström on 2016-04-04

path_main <- "/home/Crisp/clinseq/from_rasta/BREASTv6/"
#find the tumor files with deviant names
bamFiles <- dir(pattern="T*panel_v1.bam$",path=path_main,recursive=TRUE) 
bamFiles <- bamFiles[-grep("N_panel_v1",bamFiles)]
bamFiles <- bamFiles[-grep("T_panel_v1",bamFiles)]
bamFiles_Path <- paste(paste(path_main,bamFiles,sep=""),collapse=" ")
output_dir <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/batchResultsDevs/"
setwd(output_dir)
ref_name <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/batchResults/referenceFlatPanelV1.cnn"

#run the tumors 
systemCommand <- paste("cnvkit.py batch ", bamFiles_Path, " -r ", ref_name, " -p 5 --output-dir ", output_dir,sep="")
system(systemCommand)



