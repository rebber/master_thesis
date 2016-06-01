#Installera från github
remove.packages("QDNAseq")
install.packages("devtools")
library("devtools")
remove.packages("Rsamtools")
source("http://bioconductor.org/biocLite.R") #OBS! Installerar Bioconductor 3.0 pga gammal R version.
biocLite("Rsamtools") #För Rsamtools 1.22 behövs Bioconductor 3.2, därför endast Rsamtools 1.18.3 installeras?
install_github("Bioconductor-mirror/Rsamtools") #Ger error trots att gamla package tas bort först: "ERROR: compilation failed for package ‘Rsamtools’* removing ‘/home/rebecka.bergstrom/R/x86_64-redhat-linux-gnu-library/3.1/Rsamtools’"
install_github("ccagc/QDNAseq")

library("QDNAseq") #måste ladda om library om det redan finns laddat sedan tidigare, annars ladda det för första gången

instpack <- installed.packages() #lists all installed packages
grep("Rsamtools",instpack,value=TRUE)
