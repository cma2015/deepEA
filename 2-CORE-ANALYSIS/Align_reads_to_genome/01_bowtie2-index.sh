echo "genome: ${1}"

if [ ! -f "${1}.1.bt2" ]; then 
    bowtie2-build ${1} ${1}
fi 
