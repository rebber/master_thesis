#Save MAPD and mrpb matrices for wgs data for the Bv6 test samples
#Created by Rebecka Bergström on 2016-04-19

path_main <- "/home/Crisp/clinseq/from_rasta/BREASTv6/"
path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/"
load("/home/rebecka.bergstrom/cnvkit/BREASTv6_all/test50/test50")
samples <- gsub("_panel_v1","",testset[,1])
files_wgs <- dir(pattern="qdnaseq.txt$", path=path_main,recursive=TRUE)
files_wgs_test <- rep(NA, length(samples))
for (i in 1:length(samples)) {
  files_wgs_test[i] <- files_wgs[grep(paste("_",samples[i],sep=""),files_wgs)]
  #print(paste("Sample", i, "has", length(grep(paste(samples[i],"_wgs.qdnaseq.txt",sep=""),files_wgs)), "matches."))
  print(paste("Sample", i, "has", length(grep(paste("_",samples[i],sep=""),files_wgs)), "matches."))
}

#pre allocate matrices
MAPD_wgs <- rep(NA,length(samples))
mrpb_wgs <- rep(NA,length(samples))
MAPD_wgs_alt <- rep(NA,length(samples))
mrpb_wgs_alt <- rep(NA,length(samples))

#calculate MAPD & mrpb for each sample
for (i in 1:length(samples)) {
  dat_wgs <- read.table(paste(path_main,files_wgs_test[i],sep=""),header=TRUE,as.is=TRUE)
  dat_wgs$copynumber[which(dat_wgs$copynumber<0)] <- 1e-310
  mrpb_wgs[i] <- median(dat_wgs$readcount[which(dat_wgs$copynumber!="NA")])
  MAPD_wgs[i] <- median(abs(diff(log2(dat_wgs$copynumber[which(dat_wgs$copynumber!="NA")]))))
  mrpb_wgs_alt[i] <- median(dat_wgs$readcount)
  MAPD_wgs_alt[i] <- median(abs(diff(log2(dat_wgs$copynumber))),na.rm=TRUE)
}

layout(1,1,1)
plot(mrpb_wgs,MAPD_wgs)
points(mrpb_wgs_alt,MAPD_wgs_alt,col="red") 
# --> it makes more sense to use mrpb_wgs & MAPD_wgs than mrpb_wgs_alt & MAPD_wgs_alt, since the first more correctly represents the mrpb of used bins 
#(the alt alternative shows lower mrpb, thus making it seem closer to the ideal curve than it actually is)

write.table(mrpb_wgs,file=paste(path_CNVkit,"MAPDcomp/mrpb_wgs.txt",sep=""),quote=FALSE,sep="\t")
write.table(MAPD_wgs,file=paste(path_CNVkit,"MAPDcomp/MAPD_wgs.txt",sep=""),quote=FALSE,sep="\t")


