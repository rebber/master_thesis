# Check the lengths between targets
#Creataed by Rebecka Bergström on 2016-03-10

#Import chromosome sizes
chr_sizes <- read.table("/home/rebecka.bergstrom/cnvkit/hg19.chrom.sizes.txt",as.is=TRUE) #length of each chromosome (hg19)
chr_sizes$V1 <- gsub("chr","",chr_sizes$V1)
chr_sizes <- chr_sizes[-grep("_",chr_sizes$V1),]
chr_sizes <- chr_sizes[-grep("M",chr_sizes$V1),]

#get the target positions and calculate the lengths between them
targets_bed <- read.table("/home/rebecka.bergstrom/cnvkit/MSI-PILOT2/100kbATbins/full_tumors/batchResults/big_design.targets.slopped.target.bed", as.is=TRUE)
between_lengths <- targets_bed[2:length(targets_bed[,1]),2]-targets_bed[1:(length(targets_bed[,1])-1),3]
for (i in 1:length(between_lengths)) {
  if (between_lengths[i]<0) {
    chr_length <- chr_sizes[which(chr_sizes[,1]==targets_bed[i,1]),2]
    between_lengths[length(between_lengths)+1] <- chr_length - targets_bed[i,3] #the length between end of chromosome and last target, append to end of between_lengths
    between_lengths[i] <- targets_bed[i+1,2] - 1 #the length between first target and start of chromosome, exchange the negative length towards this
  }
}
between_lengths[length(between_lengths)+1] <- targets_bed[1,2] - 1 #the length between the first target and start of the first chromosome, append to end of between_lengths
between_lengths[length(between_lengths)+1] <- chr_sizes[21,2] - targets_bed[length(targets_bed[,3]),3] #the length between end of Y chromosome and last target, append to end of between_lengths

median(between_lengths)
mean(between_lengths)
quantile(between_lengths,probs=seq(0,1,0.125))
quantile(between_lengths[which(between_lengths>= 1e6)],probs=seq(0,1,0.125))
length(which(between_lengths>=1e6))
length(which(between_lengths>=1e6))/length(between_lengths)
length(which(between_lengths<0))/length(between_lengths)
plot(density(between_lengths[which(between_lengths>=0)]))
plot(density(between_lengths),xlim=c(0,4e6))
length(which(between_lengths>=0 & between_lengths<62500))
