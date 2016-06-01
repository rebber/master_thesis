#check targets in panel v1, v4 and big_design
#Created by Rebecka Bergström on 2016-05-10

#total genome length
chr_sizes <- read.table("/home/rebecka.bergstrom/cnvkit/hg19.chrom.sizes.txt",as.is=TRUE) #length of each chromosome (hg19)
chr_sizes$V1 <- gsub("chr","",chr_sizes$V1)
chr_sizes <- chr_sizes[-grep("_",chr_sizes$V1),]
chr_sizes <- chr_sizes[-grep("M",chr_sizes$V1),]
tot_genome_length <- sum(as.numeric(chr_sizes$V2))

#v1
bed_cnvkit_v1 <- read.table("~/cnvkit/BREASTv6_all/batchResults/clinseq_v1_targets.targets.slopped.target.bed",as.is=TRUE)
genes_per_targ_v1 <- bed_cnvkit_v1[,4]
num_targ_v1 <- length(genes_per_targ_v1)
num_no_gene_v1 <- length(which(genes_per_targ_v1=="-"))
num_no_gene_v1/num_targ_v1
diff_genes_v1 <- unique(genes_per_targ_v1)

bed_clinseq_v1 <- read.table("/proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/clinseq_v1_targets.targets.slopped.bed",as.is=TRUE)
lengths_v1 <- bed_clinseq_v1[,3] - bed_clinseq_v1[,2]
tot_lengths_v1 <- sum(lengths_v1)
part_v1 <- tot_lengths_v1/tot_genome_length
quantile(lengths_v1,probs=seq(0,1,0.1))
hist(lengths_v1)
plot(lengths_v1)

#v4
bed_clinseq_v4 <- read.table("/proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/clinseq_v4.targets.slopped.bed",as.is=TRUE)
lengths_v4 <- bed_clinseq_v4[,3] - bed_clinseq_v4[,2]
tot_lengths_v4 <- sum(lengths_v4)
part_v4 <- tot_lengths_v4/tot_genome_length
quantile(lengths_v4,probs=seq(0,1,0.1))
hist(lengths_v4)
plot(lengths_v4)

#big_design
bed_clinseq_bd <- read.table("/proj/b2010040/private/nobackup/autoseq-genome/intervals/targets/big_design.targets.slopped.bed",as.is=TRUE)
lengths_bd <- bed_clinseq_bd[,3] - bed_clinseq_bd[,2]
tot_lengths_bd <- sum(lengths_bd)
part_bd <- tot_lengths_bd/tot_genome_length
quantile(lengths_bd,probs=seq(0,1,0.1))
hist(lengths_bd)
plot(lengths_bd)
