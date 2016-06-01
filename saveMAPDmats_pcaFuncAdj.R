#save MAPD and mrpb mats for pca function adjusted log2 values. 
#Created by Rebecka Bergström on 2016-04-19

load("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/Rplots/smooth/randTestSet/pca/withoutY/")
#calculate the nonadjusted errors i the bins in Y chrom and add to the function adjusted matrix
nonadj_err_Ybins <- err_mat_refFlat[as.integer(testset[,2]),Ybins]
nonadj_log2_Ybins <- nonadj_err_Ybins + smooth_mat_refFlat[as.integer(testset[,2]),Ybins]
func_adj_log2_all_incl_Ybins <- cbind(func_adj_log2_all, nonadj_log2_Ybins)

MAPD_func_adj <- matrix(nrow=length(testset[,2]),ncol=2)
MAPD_func_adj[,1] <- apply(abs(t(diff(t(func_adj_log2_all_incl_Ybins[,which(dat_cnr_tumors_refFlat$targetType=="T")])))),1,median)
MAPD_func_adj[,2] <- apply(abs(t(diff(t(func_adj_log2_all_incl_Ybins[,which(dat_cnr_tumors_refFlat$targetType=="AT")])))),1,median)

#save the MAPD vaues 
write.table(MAPD_func_adj,file=paste(path_CNVkit,"MAPDcomp/MAPD_refFlat_func_adj.txt",sep=""),quote=FALSE,sep="\t")

#copy the median reads per bin table (since read count per bin is the same whatever adjustments we do or not)
systemCommand <- paste("cp ", path_CNVkit, "MAPDcomp/mrpb_refFlat_nonadj.txt", " ", path_CNVkit, "MAPDcomp/mrpb_refFlat_func_adj.txt",sep="")
system(systemCommand)


