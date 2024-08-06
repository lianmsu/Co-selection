#Customized script for processing diamond results
#!/bin/sh

diamond_res="~/meta/22MAGs_genes_diamond/diamondres3"
id_mapping_file="~/database/SARG/1.SARG_v3.2_20220917_Full_structure.txt"
outpath="~/bingene/sarg"
mkdir -p $outpath
gene=$(basename "$outpath")
find $diamond_res -name "*ARG_out" -type f -exec sh -c 'folder=$(dirname "{}"); basename=$(basename "$folder"); awk -v foldername="$basename" '\''{print $0 "\t" foldername}'\'' "{}" ' \; > ${outpath}/merged_${gene}_output_file
awk -F'\t' 'NR==FNR{a[$2]=$4; next} {print $0 "\t" (a[$2] ? a[$2] : "")}' "${id_mapping_file}" "${outpath}/merged_${gene}_output_file" > "${outpath}"/${gene}_gene_out1
python3 turn_to_table.py ${gene}_gene_out1

id_mapping_file="~/database/MRG_structure.txt"
outpath="~/bingene/mrg"
mkdir -p $outpath
gene=$(basename "$outpath")
find $diamond_res -name "*MRG_out" -type f -exec sh -c 'folder=$(dirname "{}"); basename=$(basename "$folder"); awk -v foldername="$basename" '\''{print $0 "\t" foldername}'\'' "{}" ' \; > ${outpath}/merged_${gene}_output_file
awk -F'\t' 'NR==FNR{a[$2]=$4; next} {print $0 "\t" (a[$2] ? a[$2] : "")}' "${id_mapping_file}" "${outpath}/merged_${gene}_output_file" > "${outpath}"/${gene}_gene_out1
python3 turn_to_table.py ${gene}_gene_out1

id_mapping_file="~/database/ADG_structure.txt"
outpath="~/bingene/adg"
mkdir -p $outpath
gene=$(basename "$outpath")
find $diamond_res -name "*ADG_out" -type f -exec sh -c 'folder=$(dirname "{}"); basename=$(basename "$folder"); awk -v foldername="$basename" '\''{print $0 "\t" foldername}'\'' "{}" ' \; > ${outpath}/merged_${gene}_output_file
awk -F'\t' 'NR==FNR{a[$2]=$4; next} {print $0 "\t" (a[$2] ? a[$2] : "")}' "${id_mapping_file}" "${outpath}/merged_${gene}_output_file" > "${outpath}"/${gene}_gene_out1
python3 turn_to_table.py ${gene}_gene_out1


diamondres=/home/LPW/meta/22MAGs_genes_diamond/diamondres3
cd /home/LPW/meta/22MAGs_genes_diamond/diamondres3
echo -e "genome\tARG\tMRG\tAromaDeg\tMGE" >> tongji2.txt

for forder in *
do
echo ${forder}
cd ${forder}

#得到每个文件的行数
ARG=$(awk 'END{print NR}' *ARG_out)
MRG=$(awk 'END{print NR}' *MRG_out)
AromaDeg=$(awk 'END{print NR}' *AromaDeg_out)
MGE=$(awk 'END{print NR}' *MGE_out)

echo -e "${forder}\t${ARG}\t${MRG}\t${AromaDeg}\t${MGE}" >> ${diamondres}/tongji2.txt

cd /home/LPW/meta/22MAGs_genes_diamond/diamondres3
done
