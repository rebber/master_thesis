#Run CNVkit for BREASTv6 tumors, test set, use all normals as reference
#Created by Rebecka Bergström on 2016-04-11

path_main <- "/home/Crisp/clinseq/from_rasta/BREASTv6/"
load("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/test50/test50") #loads the testset variable
bamFiles <- c() 
for (i in 1:length(testset[,1])) {
  bamFile <- dir(pattern=paste(testset[i,1],".bam$",sep=""),path=path_main,recursive=TRUE)
  bamFiles <- append(bamFiles,bamFile)
}
bamFiles_Path <- paste(paste(path_main,bamFiles,sep=""),collapse=" ")

#All normals as ref
bamFiles_ref_Path <- paste(path_main,"*/datareports/*/*N_panel_v1.bam",sep="") # tumor files
#normals <- dir(pattern="N_panel_v1.bam$",path_main,recursive=TRUE) 
#length(unique(basename(file_path_sans_ext(normals)))) #how many unique normals? 327 st. 
ref_name <- "reference327normalsV1.cnn"
output_dir <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/ref327norm/batchResults/"
setwd(output_dir)

#run the tumors 
systemCommand <- paste("cnvkit.py batch ", bamFiles_Path, " --normal ", bamFiles_ref_Path, " --targets /proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/clinseq_v1_targets.targets.slopped.bed ",
                       "--split --annotate ~/cnvkit/refFlat_reformatted.txt --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta ",
                       "--access ~/cnvkit/access-5kb-mappable.hg19_reformatted.bed --output-reference ", ref_name, " ",
                       "--output-dir ",output_dir, sep="")
system(systemCommand)



