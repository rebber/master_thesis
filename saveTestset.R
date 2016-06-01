#Randomly select 50 samples to be test set for pca analysis
#Created by Rebecka Bergström on 2016-04-08

path_CNVkit <- "/home/rebecka.bergstrom/cnvkit/BREASTv6_all/"
files_cnr <- dir(pattern="T*panel_v1.txt", paste(path_CNVkit,"cnrRcComb",sep=""), recursive = TRUE)
samples <- unique(basename(file_path_sans_ext(files_cnr)))

testset <- sample(samples,50)
testset_dev <- testset[-grep("T_p",testset)]
length(testset_dev)/length(testset)
53/length(samples)
ind <- c()
for (i in 1:length(testset)) {
  ind <- append(ind,which(testset[i]==samples))
}

o <- order(ind)
testset <- cbind(testset[o], ind[o])
save(testset, file="~/cnvkit/BREASTv6_all/test50/test50")
