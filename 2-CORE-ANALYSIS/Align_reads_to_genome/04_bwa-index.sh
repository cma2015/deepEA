echo "genome: ${1}" 

if [ ! -f "${1}.bwt" ]; then 
    /home/miniconda2/bin/bwa index ${1}
fi 