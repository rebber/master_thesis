setwd("/home/Crisp/rebecka/")


#Installera från github
remove.packages("QDNAseq")
install.packages("devtools")
library(devtools)
remove.packages("Rsamtools")
source("http://bioconductor.org/biocLite.R") 
biocLite("Rsamtools") 
library(Rsamtools)
install_github("Bioconductor-mirror/Rsamtools") #Ger error trots att gamla package tas bort först: "ERROR: compilation failed for package ‘Rsamtools’* removing ‘/home/rebecka.bergstrom/R/x86_64-redhat-linux-gnu-library/3.1/Rsamtools’"
install_github("ccagc/QDNAseq")

library("QDNAseq") #måste ladda om library om det redan finns laddat sedan tidigare, annars ladda det för första gången

instpack <- installed.packages() #lists all installed packages
grep("Rsamtools",instpack,value=TRUE)

