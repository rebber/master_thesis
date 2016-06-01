#Run CNVkit for BREASTv6 tumors, test set, use 5 normals as reference
#Created by Rebecka Bergström on 2016-04-09

path_main <- "/home/Crisp/clinseq/from_rasta/BREASTv6/"
load("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/test50/test50") #loads the testset variable
bamFiles <- c() 
for (i in 1:length(testset[,1])) {
  bamFile <- dir(pattern=paste(testset[i,1],".bam$",sep=""),path=path_main,recursive=TRUE)
  bamFiles <- append(bamFiles,bamFile)
}
bamFiles_Path <- paste(paste(path_main,bamFiles,sep=""),collapse=" ")

#5 random normals as ref
#files <- dir(pattern="N_panel_v1.bam$", path_main, recursive = TRUE)
#samples <- unique(basename(file_path_sans_ext(files)))
#normals <- sample(samples,5)
#got 5 random normals, now use them always (not 5 new each time this may be run)
normals <- c("norm1N_panel_v1", "norm2N_panel_v1", "norm3N_panel_v1",  "norm4N_panel_v1", "norm5N_panel_v1") 
#exchange "norm1" etc for actual sample IDs
bamFiles_ref <- c() 
for (i in 1:length(normals)) {
  bamFile_ref <- dir(pattern=paste(normals[i],".bam$",sep=""),path=path_main,recursive=TRUE)
  bamFiles_ref <- append(bamFiles_ref,bamFile_ref)
}
bamFiles_ref_Path <- paste(paste(path_main,bamFiles_ref,sep=""),collapse=" ")
ref_name <- "reference5normalsV1.cnn"
output_dir <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/ref5norm/batchResults/"
setwd(output_dir)

#run the tumors 
systemCommand <- paste("cnvkit.py batch ", bamFiles_Path, " --normal ", bamFiles_ref_Path, " --targets /proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/clinseq_v1_targets.targets.slopped.bed ",
                       "--split --annotate ~/cnvkit/refFlat_reformatted.txt --fasta /proj/b2010040/private/nobackup/autoseq-genome/genome/human_g1k_v37_decoy.fasta ",
                       "--access ~/cnvkit/access-5kb-mappable.hg19_reformatted.bed --output-reference ", ref_name, " ",
                       "--output-dir ",output_dir, sep="")
system(systemCommand)



