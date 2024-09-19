# 16S Usearch
chmod +x usearch

./usearch -fastq_mergepairs ./*_R1.raw.fastq -relabel @ -fastqout merged.fq

./usearch -fastx_truncate merged.fq -stripleft 27 -stripright 28 -fastqout stripped.fq 
./usearch -fastq_filter stripped.fq -fastq_maxee 1.0 -fastaout filtered.fa -relabel Filt 
./usearch -fastx_uniques filtered.fa -sizeout -relabel Uniq -fastaout uniques.fa 

./usearch -unoise3 uniques.fa -zotus asv.fa

head -n4 asv.fa 
sed 's/Zotu/ASV_/g'  asv.fa >  otusasv.fa
head -n10 otusasv.fa 
 cp -f  otusasv.fa otus.fa 
head -n10 otus.fa 

./usearch -otutab merged.fq -db otus.fa -otutabout otutab_raw.txt 

./usearch -otutab_rare otutab_raw.txt -sample_size 50000 -output otutab.txt
./usearch -makeudb_usearch SINTAX_MiDAS4.8.1.fa -output SINTAX_MiDAS4.8.1.udb 
./usearch11 -sintax otus.fa -db SINTAX_MiDAS4.8.1.udb -strand both -tabbedout sintax.txt -sintax_cutoff 0.5
 
