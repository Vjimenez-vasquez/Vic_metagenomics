# Vic_metagenomics
collection of codes for metagenomics abundance table

## Usage 
```r
#!/bin/bash
# -*- ENCODING: UTF-8 -*-.
#1#identificacion con kraken#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz
kraken2 --paired --use-names --db /home/vjimenez/Documentos/minikraken2_v2_8GB_201904_UPDATE --report-zero-counts --threads 15 $r1 $r2 --report ${prefix}_report.txt --output ${prefix}_kraken2.out ;
done ;
mkdir kraken_run2 ; 
mv *.out *.txt kraken_run2/ ;

#2#filtrado con el genoma humano#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 15 human.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 15 -bS -T human.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 15 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 15 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 15 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 15 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 15 ${prefix}.bam ;

#4# remover los archivos intermediarios#
rm ${prefix}_uno.bam ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;

#5# extraer del archivo .bam, todos los reads alineados con la referencia y generar dos fastq (f y r)#
bam2fastq --no-aligned -o ${prefix}_unal#.fastq ${prefix}.bam ;
done ;

#6# mover estos nuevos fastq a una carpeta nueva para prerar el ensamblaje de novo#
mkdir human_genome ; 
mv human.fasta *.amb *.ann *.bwt *.pac *.sa *.fai human_genome/ ; 
mkdir aligned_for_denovo ;
mv *.fastq aligned_for_denovo/ ;
cd aligned_for_denovo ;

#7# generar archivos fastq solo con reads pareados#
for r1 in *fastq
do
prefix=$(basename $r1 _unal_1.fastq)
r2=${prefix}_unal_2.fastq
fastq_pair $r1 $r2 ; 
rm *.single.fq ;
done ;

#8# realizar en ensamblaje de novo con spades#
for r1 in *fq
do
prefix=$(basename $r1 _unal_1.fastq.paired.fq)
r2=${prefix}_unal_2.fastq.paired.fq
spades --pe1-1 $r1 --pe1-2 $r2 --careful -t 15 -o ${prefix}_spades ;
mv ${prefix}_spades/scaffolds.fasta ${prefix}_spades/${prefix}_spades_scaffolds.fasta ;
mv ${prefix}_spades/${prefix}_spades_scaffolds.fasta . ;
done ;
rmdir *fq_spades ;
mv *scaffolds.fasta .. ;
cd .. ;
mkdir scaffolds ; 
mv *scaffolds.fasta script_2.R scaffolds/ ;
cd scaffolds/ ;

#9# now merge all contigs and change names according to sample#
Rscript script_2.R ;
mv *_1_spades_scaffolds.fasta script_2.R .. ;
cd .. ; 
mkdir CAT_identification ; 
cp *_1_spades_scaffolds.fasta CAT_identification/ ;
cd CAT_identification/ ;
cat *_1_spades_scaffolds.fasta > all.fasta ; 

#10# CAT/BAT#
CAT contigs -c all.fasta -d /home/vjimenez/Descargas/CAT_prepare_20210107/2021-01-07_CAT_database -t /home/vjimenez/Descargas/CAT_prepare_20210107/2021-01-07_taxonomy -o all2 ;
CAT add_names -i all2.contig2classification.txt -o all2.contig2classification.official_names.txt -t /home/vjimenez/Descargas/CAT_prepare_20210107/2021-01-07_taxonomy --only_official ;
CAT summarise -c all.fasta -i all2.contig2classification.official_names.txt -o all2.summary.txt ;
cd .. ;

#11# looping indexing#
for t1 in *scaffolds.fasta
do
prefix=$(basename $t1 _1_spades_scaffolds.fasta)
bwa index $t1 ;
done ;

#12# looping mapping#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz
r3=${prefix}_1_spades_scaffolds.fasta

#13# instrucciones para generar el archivo .bam#
bwa mem -t 15 $r3 $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 15 -bS -T $r3 ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 15 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 15 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 15 ${prefix}_tresa.bam -o ${prefix}_count.bam ;
samtools index -@ 15 ${prefix}_count.bam ;

#14# remover los archivos intermediarios#
rm ${prefix}_uno.bam ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ;
done ;

#15# loop estimates read counts per contig#
for p1 in *_count.bam
do
prefix=$(basename $p1 _count.bam)
samtools view -@ 15 $p1 | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c > ${prefix}_readcounts.txt ;
done ;
rm *.amb *.ann *.bwt *.fai *.pac *.sa ;
cat *_readcounts.txt > counts.txt ;
grep "NODE" counts.txt > counts2.txt ;

#16# merge odentification and abundances#
Rscript script_3.R ; 
exit
``` 
