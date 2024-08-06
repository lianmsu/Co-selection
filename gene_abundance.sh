#!/bin/sh

source activate args_oap
ipt=~/22ARG-JILI/in

DB=~/22ARG-JILI/database2/4.SARG_v3.2_20220917_Short_subdatabase.fasta
struc=~/22ARG-JILI/database2/4.SARG_v3.2_20220917_Short_subdatabase_structure.txt
mkdir ~/22ARG-JILI/argjiliout2
echo Step1
args_oap stage_one -i ${ipt}  -o ~/22ARG-JILI/argjiliout2  -t 28 --database ${DB}
echo Step2
args_oap stage_two -i ~/22ARG-JILI/argjiliout2 --database ${DB} -t 28 --structure1 ${struc}

DB=~/meta/database/AromaDeg/AromaDeg.fasta
struc=~/meta/aromade/AromaDeg2_structure.txt
mkdir ~/ADG
echo Step1
args_oap stage_one -i ${ipt}  -o ~/ADG  -t 28 --database ${DB}
echo Step2
args_oap stage_two -i ~/ADG --database ${DB} -t 28 --structure1 ${struc}

DB=~/meta/database/MRG/MRG.fasta
struc=~/meta/database/MRG/MRG_structure.txt
mkdir ~/MRG
echo Step1
args_oap stage_one -i ${ipt}  -o ~/MRG  -t 28 --database ${DB}
echo Step2
args_oap stage_two -i ~/MRG --database ${DB} -t 28 --structure1 ${struc}

DB=~/meta/database/MGE/MGE.fasta
struc=~/meta/database/MGE/MGE_structure.txt
mkdir ~/MGE
echo Step1
args_oap stage_one -i ${ipt}  -o ~/MGE  -t 28 --database ${DB}
echo Step2
args_oap stage_two -i ~/MGE --database ${DB} -t 28 --structure1 ${struc}





