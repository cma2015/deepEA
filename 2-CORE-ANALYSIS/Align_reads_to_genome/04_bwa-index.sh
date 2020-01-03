echo "genome: ${1}" 

if [ ! -f "${1}.bwt" ]; then 
    bwa index ${1}
fi 