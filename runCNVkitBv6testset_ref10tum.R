#Run CNVkit for BREASTv6 tumors, test set, use 10 flat tumors as reference
#Created by Rebecka Bergström on 2016-04-09

path_main <- "/home/Crisp/clinseq/from_rasta/BREASTv6/"
load("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/test50/test50") #loads the testset variable
bamFiles <- c() 
for (i in 1:length(testset[,1])) {
  bamFile <- dir(pattern=paste(testset[i,1],".bam$",sep=""),path=path_main,recursive=TRUE)
  bamFiles <- append(bamFiles,bamFile)
}
bamFiles_Path <- paste(paste(path_main,bamFiles,sep=""),collapse=" ")
bamFiles_ref <- dir(pattern="T.*panel_v1.bam","/home/Crisp/rebecka/breastv7_bams_with_no_cna/",recursive=TRUE) #11 flat tumors
bamFiles_ref <- bamFiles_ref[-grep("a1",bamFiles_ref)] #remove the least flat one (exchange a1 towards actual sample ID)
bamFiles_ref_Path <- paste("/home/Crisp/rebecka/breastv7_bams_with_no_cna/",bamFiles_ref," ",sep="",collapse = "")
ref_name <- "reference10tumorsV1.cnn"
output_dir <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/ref10tum/batchResults/"
setwd(output_dir)

#run the tumors 
systemCommand <- paste("cnvkit.py batch ", bamFiles_Path, " --normal ", bamFiles_ref_Path, "--targets /proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/clinseq_v1_targets.targets.slopped.bed ",
                       "--split --annotate ~/cnvkit/refFlat_reformatted.txt --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta ",
                       "--access ~/cnvkit/access-5kb-mappable.hg19_reformatted.bed --output-reference ", ref_name, " ",
                       "--output-dir ",output_dir, sep="")
system(systemCommand)



