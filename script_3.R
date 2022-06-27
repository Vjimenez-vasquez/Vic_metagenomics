#script to modify the counts output#

table <- read.csv("counts2.txt", header=FALSE) ;
names(table) <- "cabeza" ;
head(table) ;

b <- c()
c <- c()
e <- c()
f <- c()
#need to create a specific list names for your samples#
list <- c("Adeno-Inv","Adeno-BIO")
for (i in 1:length(list)){
b <- grep(list[i],table$cabeza, value = TRUE)
c <- gsub(paste0(" ",list[i]),paste0(",",list[i]),b)
e <- gsub(" ","",c)
f <- append(f,e)
}
g <- as.data.frame(f)
class(g)
dim(g)

library(tidyr)
h <- separate(g,"f",into=c("n_reads","contig"),sep=",")
write.csv(h,"counts3.csv", row.names=FALSE)

#made in R#
c <- read.csv("CAT_identification/all2.contig2classification.official_names.txt", header=TRUE, sep="\t") ;
names(c) <- c("contig","classification","reason","lineage","lineage.scores","superkingdom","phylum","class","order","family","genus","species")
n <- read.csv("counts3.csv", header=FALSE) ;
names(n) <- c("n_reads","contig") ;
m <- merge(c,n,by="contig", all.x=TRUE) ;
write.csv(m,"CAT_identified_abundances.csv", row.names=FALSE) ;
q("no") ;
