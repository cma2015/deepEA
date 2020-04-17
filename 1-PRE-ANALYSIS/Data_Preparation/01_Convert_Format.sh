#!/bin/bash
echo "Script:$0";
echo "SRA ID:$1";
echo "minReadLen:$2";



dirpath=/home/DeepEA/galaxy/tools/1-PRE-ANALYSIS/Data_Preparation/$3/

/home/miniconda2/bin/fasterq-dump $1 -M $2 --split-3 -O $dirpath ;


count=$(ll ${dirpath} | wc -l)
echo $count
#ls -l ${dirpath} | wc -l > $6



name=$3
if [ ${count} == "2" ] ; then
    mv ${dirpath}*_1.fastq $4 ; #cat ${dirpath}*_1.fastq > $4 ;
    mv ${dirpath}*_2.fastq $5 ;
else
    mv ${dirpath}*.fastq  $4 ;
fi

