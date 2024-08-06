#!/bin/sh

orf=~/meta/21MAGs_prodigal/MAGs_gene
out=~/meta/22MAGs_genes_diamond/diamondres3
#database
ARGdb=~/meta/database/SARG/SARG_v3.2_diamond.dmnd
MRGdb=~/meta/database/MRG/MRG.dmnd
MGEdb=~/meta/database/MGE/MGE.dmnd
AromaDegdb=~/meta/database/AromaDeg/AromaDeg.fasta.dmnd

cd ${orf}
for folder in *
do
mkdir -p ${out}/${folder}
query=${folder}/protein_seq.fa
diamond blastp -k 1 -e 0.00001 -p 14 -d ${ARGdb} -q ${query} -o ${out}/${folder}/${folder}_ARG_out --id 80 --query-cover 80
diamond blastp -k 1 -e 0.00001 -p 14 -d ${MGEdb} -q ${query} -o ${out}/${folder}/${folder}_MGE_out --id 80 --query-cover 80
diamond blastp -k 1 -e 0.00001 -p 14 -d ${MRGdb} -q ${query} -o ${out}/${folder}/${folder}_MRG_out --id 80 --query-cover 80
diamond blastp -k 1 -e 0.00001 -p 14 -d ${AromaDegdb} -q ${query} -o ${out}/${folder}/${folder}_AromaDeg_out --id 80 --query-cover 80
done
