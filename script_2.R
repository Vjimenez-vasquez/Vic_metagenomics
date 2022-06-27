#9# now change fasta headers and merge all data#
library(seqinr) ;
r <- dir() ;
head <- gsub("_spades.*","",r) ;
a <- 0 ;
for (i in 1:length(head)){ ;
a <- read.fasta(r[i]) ;
names(a) <- paste0(rep(head[i],length(a)),"_",names(a)) ;
write.fasta(a,names(a), file.out=paste0(head[i],"_1_spades_scaffolds.fasta")) ;
} ;
q("no") ;
