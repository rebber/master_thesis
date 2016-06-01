#Calculate MAPD and mrpb for HOPKINS1000G test samples with two diff refs:  50 normals and 50 normals+PCA function adj respectively
#Created by Rebecka Bergström on 2016-05-31

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/HOPKINS/hopkins_1000G/"

#for the 50 normals reference:
load(paste(path_CNVkit,"test30",sep="")) #load testset 
samples <- testset[,1]

files_cnr_50norm <- paste("cnrRcComb/addedSmooth/",samples,".txt",sep="")

#preallocate MAPD & mrpb matrices
MAPD_50norm <- matrix(nrow=length(samples),ncol=2)
colnames(MAPD_50norm) <- c("targets","antitargets")
rownames(MAPD_50norm) <- samples
mrpb_50norm <- matrix(nrow=length(samples),ncol=2)
colnames(mrpb_50norm) <- c("targets","antitargets")
rownames(mrpb_50norm) <- samples

for (i in 1:length(samples)) {
  dat_cnr_50norm <- read.table(paste(path_CNVkit,files_cnr_50norm[i],sep=""),header=TRUE,as.is=TRUE)
  
  #find median read depths per bin (mrpb) and store
  mrpb_50norm[i,1] <- median(dat_cnr_50norm$reads[which(dat_cnr_50norm$targetType=="T")])
  mrpb_50norm[i,2] <- median(dat_cnr_50norm$reads[which(dat_cnr_50norm$targetType=="AT")])
  
  #find the median absolute pairwise difference (MAPD) for the log2CN-ratios and store
  MAPD_50norm[i,1] <- median(abs(diff(dat_cnr_50norm$log2[which(dat_cnr_50norm$targetType=="T")])))
  MAPD_50norm[i,2] <- median(abs(diff(dat_cnr_50norm$log2[which(dat_cnr_50norm$targetType=="AT")])))
  
}

write.table(mrpb_50norm,file=paste(path_CNVkit,"MAPDcomp/mrpb_ref50norm.txt",sep=""),quote=FALSE,sep="\t")
write.table(MAPD_50norm,file=paste(path_CNVkit,"MAPDcomp/MAPD_ref50norm.txt",sep=""),quote=FALSE,sep="\t")

#################################

#for 50 norm ref + fitted function PCA adjustment of bins
load("/home/rebecka.bergstrom/cnvkit/HOPKINS/hopkins_1000G/pca/.RData")

#calculate the nonadjusted errors i the bins in Y chrom and add to the function adjusted matrix
#not necessary, if Ybins are not removed for HOPKINS1000G
#nonadj_err_Ybins <- err_mat_refFlat[as.integer(testset[,2]),Ybins]
#nonadj_log2_Ybins <- nonadj_err_Ybins + smooth_mat_refFlat[as.integer(testset[,2]),Ybins]
#func_adj_log2_all_incl_Ybins <- cbind(func_adj_log2_all, nonadj_log2_Ybins)
func_adj_log2_all_incl_Ybins <- func_adj_log2_all

MAPD_func_adj <- matrix(nrow=length(testset[,2]),ncol=2)
MAPD_func_adj[,1] <- apply(abs(t(diff(t(func_adj_log2_all_incl_Ybins[,which(dat_cnr_tumors$targetType=="T")])))),1,median)
MAPD_func_adj[,2] <- apply(abs(t(diff(t(func_adj_log2_all_incl_Ybins[,which(dat_cnr_tumors$targetType=="AT")])))),1,median)

#save the MAPD vaues 
write.table(MAPD_func_adj,file=paste(path_CNVkit,"MAPDcomp/MAPD_ref50norm_func_adj.txt",sep=""),quote=FALSE,sep="\t")

#copy the median reads per bin table (since read count per bin is the same whatever adjustments we do or not)
systemCommand <- paste("cp ", path_CNVkit, "MAPDcomp/mrpb_ref50norm.txt", " ", path_CNVkit, "MAPDcomp/mrpb_ref50norm_func_adj.txt",sep="")
system(systemCommand)


