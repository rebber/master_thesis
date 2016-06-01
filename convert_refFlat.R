refFlat <- read.table("/home/rebecka.bergstrom/cnvkit/refFlat.txt",as.is=TRUE)
refFlat_nochr <- refFlat
refFlat_nochr$V3 <- gsub("chr","",refFlat_nochr$V3)
length(grep("_",refFlat_nochr$V3)) #how many has not only number for chromosome but something like "6_qbl_hap6", in this case 2829
refFlat_nochr <- refFlat_nochr[-grep("_",refFlat_nochr$V3),]
length(grep("M",refFlat$V3)) #how many mitochondrial chr:s are there, in this case 0
write.table(refFlat_nochr,"/home/rebecka.bergstrom/cnvkit/refFlat_reformatted.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

#Convert access-5kb-mappable.hg19.bed in the same way
access <- read.table("/home/rebecka.bergstrom/cnvkit/cnvkit-master/data/access-5k-mappable.hg19.bed",as.is=TRUE)
head(access)
access_nochr <- access
access_nochr$V1 <- gsub("chr", "", access_nochr$V1)
access_nochr$V1 <- gsub("M", "MT", access_nochr$V1)
head(access_nochr,20)
sum(which(access_nochr$V1=="MT"))
length(grep("_",access_nochr$V1))
write.table(access_nochr,"/home/rebecka.bergstrom/cnvkit/access-5kb-mappable.hg19_reformatted.bed",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
