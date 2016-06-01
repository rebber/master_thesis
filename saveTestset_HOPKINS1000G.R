#Randomly select 30 samples to be test set for pca analysis for HOPKINS 1000G
#Created by Rebecka Bergström on 2016-05-30

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/HOPKINS/hopkins_1000G/"
files_cnr <- dir(pattern="HOP.*.txt", paste(path_CNVkit,"cnrRcComb/addedSmooth/",sep=""), recursive = TRUE)
samples <- unique(basename(file_path_sans_ext(files_cnr)))

testset <- sample(samples,30)
testset_norm <- testset[grep("N",testset)]
samples_norm <- samples[grep("N",samples)]
length(testset_norm)/length(testset)
length(samples_norm)/length(samples)
ind <- c()
for (i in 1:length(testset)) {
  ind <- append(ind,which(testset[i]==samples))
}

o <- order(ind)
testset <- cbind(testset[o], ind[o])
save(testset, file=paste(path_CNVkit,"test30",sep=""))
