#Analyze correlations between pca clusters and prep-data
#Created by Rebecka Bergström on 2016-04-14

#load the PCA workspace created with plotPCA_Bv6AllT.R
load("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/Rplots/smooth/randTestSet/pca/withoutY/.RData")
setwd(paste(path_CNVkit,"Rplots/smooth/randTestSet/pca/withoutY/prep_cor/",sep=""))

#########################

#organize the prep data

#read in the prep data
prep_data <- read.table("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/clinseq_breast_prep_info.csv",header=TRUE,as.is=TRUE,sep=";",dec=",")
prep_data <- prep_data[,-(25:length(prep_data[1,]))] #columns only containing NA
prep_data <- prep_data[-366,] #empty row
#correct some columns
prep_data$sample_name <- gsub("-T","T",prep_data$sample_name)
for (i in 1:length(prep_data$sample_name)) {
  if (prep_data$cap_pool_name[i]=="") {
    prep_data$cap_pool_name[i] <- prep_data$cap_pool_name[i-1]
  }
}
prep_data$cap_pool_name[(length(prep_data$cap_pool_name)-5):length(prep_data$cap_pool_name)] <- ""

#take only numeric columns
prep_data_numeric <- prep_data[,c(2,3,6:9,11,14:17,19)]

#which rows in prep_data corresponds to my training set samples
samples_train <- gsub("_panel_v1","",samples[-as.integer(testset[,2])])
train_prep_ind <- c()
for (i in 1:length(samples_train)) {
  ind <- which(samples_train[i]==prep_data$sample_name)
  if (length(ind)==1) {
    train_prep_ind <- append(train_prep_ind,ind)
  }
}

#which rows in prep_data corresponds to my test set samples
samples_test_noend <- gsub("_panel_v1","",samples_test)
test_prep_ind <- c()
for (i in 1:length(samples_test_noend)) {
  ind <- which(samples_test_noend[i]==prep_data$sample_name)
  if (length(ind)==1) {
    test_prep_ind <- append(test_prep_ind,ind)
  }
}


################################

#check correlation with numeric columns, use spearman correlation
cor_num <- cor(x_inter,prep_data_numeric[train_prep_ind,],method="s")
cor_num_pco <- cor(x_inter,prep_data_numeric[train_prep_ind,],use = "p",method="s")
length(which(is.na(prep_data_numeric$prep_vol[train_prep_ind])))
length(which(is.na(prep_data_numeric$prep_start_date[train_prep_ind])))
quantile(cor_num_pco)

############

#Show all correlations with numeric data in barplot
filename_cor_prepNum <- "cor_prepNum_PC1-4.jpeg"
jpeg(filename_cor_prepNum,width=9,height=6,quality=90,units="in",res=600)
layout(matrix(1,1,1))
par(mar=c(10,4,3,1))
barplot(cor_num_pco, las=2, beside=TRUE,ylab="Spearman correlation",names.arg=gsub("_"," ",colnames(cor_num_pco)))
legend("topleft",c(paste("PC",1:4,sep="")),fill=grey.colors(4))
dev.off()

#Show correlations with date and prep vol (which give the largest correlation coeffs)
filename_cor_prepNum_large <- "cor_prepNum_large_PC1-4.jpeg"
jpeg(filename_cor_prepNum_large,width=9,height=4,quality=90,units="in",res=600)
layout(matrix(1,1,1))
par(mar=c(5,4,3,1))
barplot(cor_num_pco[,c(7,11)], beside=TRUE,ylab="Spearman correlation",names.arg=gsub("_"," ",colnames(cor_num_pco[,c(7,11)])))
legend("topleft",c(paste("PC",1:4,sep="")),fill=grey.colors(4))
dev.off()

#####################

#Show how PC1-4 relate to prep start date 
filename_cor_prepDate <- "cor_prepDate_PC1-4.jpeg"
jpeg(filename_cor_prepDate,width=9,height=9,quality=90,units="in",res=600)
layout(matrix(1:4,2,2,byrow=TRUE))
par(mar=c(5,4.3,3.2,1.7)+0.1)
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,1],type="n",xlab="library prep date",ylab="PC1",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[1,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,1],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,1],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,1],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,1],col="black")
legend("topright",c("cluster 1", "cluster 2", "cluster 3", "outliers"),pch=1,col=c("red","green","blue","black"))
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,2],type="n",xlab="library prep date",ylab="PC2",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[2,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,2],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,2],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,2],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,2],col="black")
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,3],type="n",xlab="library prep date",ylab="PC3",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[3,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,3],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,3],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,3],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,3],col="black")
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,4],type="n",xlab="library prep date",ylab="PC4",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[4,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,4],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,4],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,4],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,4],col="black")
dev.off()



length(which(is.na(prep_data$prep_start_date[train_prep_ind])))

################

#Show how PC1-4 relate to prep vol
filename_cor_prepVol <- "cor_prepVol_PC1-4.jpeg"
jpeg(filename_cor_prepVol,width=9,height=9,quality=90,units="in",res=600)
layout(matrix(1:4,2,2,byrow=TRUE))
par(mar=c(5,4.3,3.2,1.7)+0.1)
plot(prep_data$prep_vol[train_prep_ind],x_inter[,1],type="n",xlab="out vol (µl) library prep",ylab="PC1",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[1,11],2)))
points(prep_data$prep_vol[train_prep_ind][clInd1],x_inter[clInd1,1],col="red")
points(prep_data$prep_vol[train_prep_ind][clInd2],x_inter[clInd2,1],col="green")
points(prep_data$prep_vol[train_prep_ind][clInd3],x_inter[clInd3,1],col="blue")
points(prep_data$prep_vol[train_prep_ind][clInd0],x_inter[clInd0,1],col="black")
legend("topleft",c("cluster 1", "cluster 2", "cluster 3", "outliers"),pch=1,col=c("red","green","blue","black"))
plot(prep_data$prep_vol[train_prep_ind],x_inter[,2],type="n",xlab="out vol (µl) library prep",ylab="PC2",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[2,11],2)))
points(prep_data$prep_vol[train_prep_ind][clInd1],x_inter[clInd1,2],col="red")
points(prep_data$prep_vol[train_prep_ind][clInd2],x_inter[clInd2,2],col="green")
points(prep_data$prep_vol[train_prep_ind][clInd3],x_inter[clInd3,2],col="blue")
points(prep_data$prep_vol[train_prep_ind][clInd0],x_inter[clInd0,2],col="black")
plot(prep_data$prep_vol[train_prep_ind],x_inter[,3],type="n",xlab="out vol (µl) library prep",ylab="PC3",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[3,11],2)))
points(prep_data$prep_vol[train_prep_ind][clInd1],x_inter[clInd1,3],col="red")
points(prep_data$prep_vol[train_prep_ind][clInd2],x_inter[clInd2,3],col="green")
points(prep_data$prep_vol[train_prep_ind][clInd3],x_inter[clInd3,3],col="blue")
points(prep_data$prep_vol[train_prep_ind][clInd0],x_inter[clInd0,3],col="black")
plot(prep_data$prep_vol[train_prep_ind],x_inter[,4],type="n",xlab="out vol (µl) library prep",ylab="PC4",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[4,11],2)))
points(prep_data$prep_vol[train_prep_ind][clInd1],x_inter[clInd1,4],col="red")
points(prep_data$prep_vol[train_prep_ind][clInd2],x_inter[clInd2,4],col="green")
points(prep_data$prep_vol[train_prep_ind][clInd3],x_inter[clInd3,4],col="blue")
points(prep_data$prep_vol[train_prep_ind][clInd0],x_inter[clInd0,4],col="black")
dev.off()

t.test(x_inter[which(prep_data$prep_vol[train_prep_ind]==19),2],x_inter[which(prep_data$prep_vol[train_prep_ind]==26),2])
length(which(is.na(prep_data$prep_vol[train_prep_ind])))

##################

#Check if cluster 3 samples have something in common
View(prep_data[train_prep_ind[clInd3],])
View(prep_data)
#They actually are closely positioned on all plates used and therefore the low on-target rate 
  # is likely due to some mishap in that area


##################

save.image("~/.RData")
plot3d(prep_data$prep_vol[train_prep_ind],as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,2],type="n")
points3d(prep_data$prep_vol[train_prep_ind][clInd1],as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,2],col="red")
points3d(prep_data$prep_vol[train_prep_ind][clInd2],as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,2],col="green")
points3d(prep_data$prep_vol[train_prep_ind][clInd3],as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,2],col="blue")
points3d(prep_data$prep_vol[train_prep_ind][clInd0],as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,2],col="black")


#cor between PC2 and prep_conc_after_pcr=-0.11
length(which(is.na(prep_data$prep_conc_after_pcr[train_prep_ind])))
plot(prep_data$prep_conc_after_pcr[train_prep_ind],x_inter[,2],type="n")
points(prep_data$prep_conc_after_pcr[train_prep_ind][clInd1],x_inter[clInd1,2],col="red")
points(prep_data$prep_conc_after_pcr[train_prep_ind][clInd2],x_inter[clInd2,2],col="green")
points(prep_data$prep_conc_after_pcr[train_prep_ind][clInd3],x_inter[clInd3,2],col="blue")
points(prep_data$prep_conc_after_pcr[train_prep_ind][clInd0],x_inter[clInd0,2],col="black")
legend("topleft",c("cluster 1", "cluster 2", "cluster 3", "outliers"),pch=1,col=c("red","green","blue","black"))
t.test(x_inter[which(prep_data$prep_conc_after_pcr[train_prep_ind]==19),2],x_inter[which(prep_data$prep_conc_after_pcr[train_prep_ind]==26),2])

####################

#Show how prep date of test samples and their cluster agree with training samples
layout(matrix(1:4,2,2,byrow=TRUE))
par(mar=c(5,4.3,3.2,1.7)+0.1)
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,1],type="n",xlab="library prep date",ylab="PC1",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[1,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,1],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,1],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,1],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,1],col="black")
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==1)],pc_test_all[which(assign_clust_all==1),1],col="red",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==2)],pc_test_all[which(assign_clust_all==2),1],col="green",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==3)],pc_test_all[which(assign_clust_all==3),1],col="blue",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==0)],pc_test_all[which(assign_clust_all==0),1],col="black",pch=8)
legend("topright",c("cluster 1", "cluster 2", "cluster 3", "outliers"),pch=1,col=c("red","green","blue","black"))
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,2],type="n",xlab="library prep date",ylab="PC2",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[2,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,2],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,2],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,2],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,2],col="black")
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==1)],pc_test_all[which(assign_clust_all==1),2],col="red",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==2)],pc_test_all[which(assign_clust_all==2),2],col="green",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==3)],pc_test_all[which(assign_clust_all==3),2],col="blue",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==0)],pc_test_all[which(assign_clust_all==0),2],col="black",pch=8)
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,3],type="n",xlab="library prep date",ylab="PC3",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[3,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,3],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,3],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,3],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,3],col="black")
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==1)],pc_test_all[which(assign_clust_all==1),3],col="red",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==2)],pc_test_all[which(assign_clust_all==2),3],col="green",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==3)],pc_test_all[which(assign_clust_all==3),3],col="blue",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==0)],pc_test_all[which(assign_clust_all==0=,3],col="black",pch=8)
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,4],type="n",xlab="library prep date",ylab="PC4",cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[4,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,4],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,4],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,4],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,4],col="black")
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==1)],pc_test_all[which(assign_clust_all==1),4],col="red",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==2)],pc_test_all[which(assign_clust_all==2),4],col="green",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==3)],pc_test_all[which(assign_clust_all==3),4],col="blue",pch=8)
points(as.Date(as.character(prep_data$prep_start_date[test_prep_ind]),format="%Y%m%d")[which(assign_clust_all==0)],pc_test_all[which(assign_clust_all==0),4],col="black",pch=8)


########################################

#Make combined plots, correlations with prep data and MAPD etc in the same plot.

#Plot combined barplot including correlations with MAPD etc.
cor_tot_reads <- cor(x_inter,tot_reads[-as.integer(testset[,2]),1],method="s")
cor_ontarget_rate <- cor(x_inter,tot_reads[-as.integer(testset[,2]),2],method="s")
cor_mrpb_targ <- cor(pca_refFlat$x[,1:4],as.numeric(mrpb_targets_refFlat[-as.integer(testset[,2]),2]),method="s")
cor_mrpb_antitarg <- cor(pca_refFlat$x[,1:4],as.numeric(mrpb_antitargets_refFlat[-as.integer(testset[,2]),2]),method="s")
cor_MAPD_targ <- cor(pca_refFlat$x[,1:4],as.numeric(MAPD_targets_refFlat[-as.integer(testset[,2]),2]),method="s")
cor_MAPD_antitarg <- cor(pca_refFlat$x[,1:4],as.numeric(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2]),method="s")
all_cor <- cbind(cor_MAPD_targ,cor_MAPD_antitarg,cor_mrpb_targ,cor_mrpb_antitarg,cor_tot_reads,cor_ontarget_rate,cor_num_pco[,c(7,11)])
colnames(all_cor) <- c("MAPD targets","MAPD antitargets","mrpb targets","mrpb antitargets","total #reads","on-target rate","library prep date", "library prep out vol")
filename_cor_combined <- "cor_combined_PC1-4.jpeg"
jpeg(filename_cor_combined,width=9,height=6,quality=90,units="in",res=600)
layout(matrix(1,1,1))
par(mar=c(9,4,3,1)+0.1)
barplot(all_cor,beside=TRUE,las=2,ylab="Spearman correlation")
#legend(x=14,y=0.77,paste("PC",1:4,sep=""),fill=grey.colors(4))
legend("topleft",paste("PC",1:4,sep=""),fill=grey.colors(4))
dev.off()

#Show chosen correlations for PCs through scatterplots
filename_cor_scatter_chosen <- "cor_scatter_chosen.jpeg"
jpeg(filename_cor_scatter_chosen,width=9,height=6,quality=90,units="in",res=600)
layout(matrix(1:6,nrow=2,ncol=3,byrow=TRUE))
par(mar=c(4.5,4.3,3.8,1.7)+0.1)
plot(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2],pca_refFlat$x[,1],ylab="PC1",xlab="MAPD antitargets",main=paste("correlation =",signif(cor_MAPD_antitarg[1],2)),type="n",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd1],pca_refFlat$x[clInd1,1],col="red")
points(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd2],pca_refFlat$x[clInd2,1],col="green")
points(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd3],pca_refFlat$x[clInd3,1],col="blue")
points(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd0],pca_refFlat$x[clInd0,1],col="black")
legend("topleft",c("cluster 1","cluster 2","cluster 3","outliers"),col=c("red","green","blue","black"),pch=1)
plot(tot_reads[-as.integer(testset[,2]),2],x_inter[,1],type="n",xlab="on-target rate",ylab="PC1",cex.main=1.5,cex.lab=1.5,cex.axis=1.5,main=paste("correlation =",signif(cor_ontarget_rate[1],2)))
points(tot_reads[-as.integer(testset[,2]),2][clInd1],x_inter[,1][clInd1],col="red")
points(tot_reads[-as.integer(testset[,2]),2][clInd2],x_inter[,1][clInd2],col="green")
points(tot_reads[-as.integer(testset[,2]),2][clInd3],x_inter[,1][clInd3],col="blue")
points(tot_reads[-as.integer(testset[,2]),2][clInd0],x_inter[,1][clInd0],col="black")
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,1],type="n",xlab="library prep date",ylab="PC1",cex.main=1.5,cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[1,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,1],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,1],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,1],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,1],col="black")
plot(MAPD_targets_refFlat[-as.integer(testset[,2]),2],pca_refFlat$x[,2],ylab="PC2",xlab="MAPD targets",main=paste("correlation =",signif(cor_MAPD_targ[2],2)),type="n",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd1],pca_refFlat$x[clInd1,2],col="red")
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd2],pca_refFlat$x[clInd2,2],col="green")
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd3],pca_refFlat$x[clInd3,2],col="blue")
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd0],pca_refFlat$x[clInd0,2],col="black")
plot(tot_reads[-as.integer(testset[,2]),1],x_inter[,2],type="n",xlab="total #reads",ylab="PC2",cex.main=1.5,cex.lab=1.5,cex.axis=1.5,main=paste("correlation =",signif(cor_tot_reads[2],2)))
points(tot_reads[-as.integer(testset[,2]),1][clInd1],x_inter[,2][clInd1],col="red")
points(tot_reads[-as.integer(testset[,2]),1][clInd2],x_inter[,2][clInd2],col="green")
points(tot_reads[-as.integer(testset[,2]),1][clInd3],x_inter[,2][clInd3],col="blue")
points(tot_reads[-as.integer(testset[,2]),1][clInd0],x_inter[,2][clInd0],col="black")
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,2],type="n",xlab="library prep date",ylab="PC2",cex.main=1.5,cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[2,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,2],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,2],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,2],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,2],col="black")
dev.off()

#Show chosen correlations for PCs through scatterplots and MAPD targets vs MAPD antitargets
layout(matrix(c(1,2,3,7,4,5,6,7),nrow=2,ncol=4,byrow=TRUE))
par(mar=c(4.5,4.3,3.8,1.7)+0.1)
plot(MAPD_targets_refFlat[-as.integer(testset[,2]),2],pca_refFlat$x[,2],ylab="PC2",xlab="MAPD targets",main=paste("correlation =",signif(cor_MAPD_targ[2],2)),type="n",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd1],pca_refFlat$x[clInd1,2],col="red")
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd2],pca_refFlat$x[clInd2,2],col="green")
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd3],pca_refFlat$x[clInd3,2],col="blue")
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd0],pca_refFlat$x[clInd0,2],col="black")
legend("bottomleft",c("cluster 1","cluster 2","cluster 3","outliers"),col=c("red","green","blue","black"),pch=1)
plot(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2],pca_refFlat$x[,1],ylab="PC1",xlab="MAPD antitargets",main=paste("correlation =",signif(cor_MAPD_antitarg[1],2)),type="n",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd1],pca_refFlat$x[clInd1,1],col="red")
points(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd2],pca_refFlat$x[clInd2,1],col="green")
points(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd3],pca_refFlat$x[clInd3,1],col="blue")
points(MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd0],pca_refFlat$x[clInd0,1],col="black")
plot(tot_reads[-as.integer(testset[,2]),1],x_inter[,2],type="n",xlab="total #reads",ylab="PC2",cex.main=1.5,cex.lab=1.5,cex.axis=1.5,main=paste("correlation =",signif(cor_tot_reads[2],2)))
points(tot_reads[-as.integer(testset[,2]),1][clInd1],x_inter[,2][clInd1],col="red")
points(tot_reads[-as.integer(testset[,2]),1][clInd2],x_inter[,2][clInd2],col="green")
points(tot_reads[-as.integer(testset[,2]),1][clInd3],x_inter[,2][clInd3],col="blue")
points(tot_reads[-as.integer(testset[,2]),1][clInd0],x_inter[,2][clInd0],col="black")
plot(tot_reads[-as.integer(testset[,2]),2],x_inter[,1],type="n",xlab="on-target rate",ylab="PC1",cex.main=1.5,cex.lab=1.5,cex.axis=1.5,main=paste("correlation =",signif(cor_ontarget_rate[1],2)))
points(tot_reads[-as.integer(testset[,2]),2][clInd1],x_inter[,1][clInd1],col="red")
points(tot_reads[-as.integer(testset[,2]),2][clInd2],x_inter[,1][clInd2],col="green")
points(tot_reads[-as.integer(testset[,2]),2][clInd3],x_inter[,1][clInd3],col="blue")
points(tot_reads[-as.integer(testset[,2]),2][clInd0],x_inter[,1][clInd0],col="black")
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,1],type="n",xlab="library prep date",ylab="PC1",cex.main=1.5,cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[1,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,1],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,1],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,1],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,1],col="black")
plot(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d"),x_inter[,2],type="n",xlab="library prep date",ylab="PC2",cex.main=1.5,cex.axis=1.5,cex.lab=1.5,main=paste("correlation =",signif(cor_num_pco[2,7],2)))
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd1],x_inter[clInd1,2],col="red")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd2],x_inter[clInd2,2],col="green")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd3],x_inter[clInd3,2],col="blue")
points(as.Date(as.character(prep_data$prep_start_date[train_prep_ind]),format="%Y%m%d")[clInd0],x_inter[clInd0,2],col="black")
plot(MAPD_targets_refFlat[-as.integer(testset[,2]),2],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2],type="n",xlab="MAPD targets",ylab="MAPD antitargets",cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd1],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd1],col="red")
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd2],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd2],col="green")
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd3],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd3],col="blue")
points(MAPD_targets_refFlat[-as.integer(testset[,2]),2][clInd0],MAPD_antitargets_refFlat[-as.integer(testset[,2]),2][clInd0],col="black")



###############











