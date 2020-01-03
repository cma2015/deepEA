echo "genome: ${1}"

if [ ! -f "${1}.1.ht2" ]; then
    /home/miniconda2/bin/hisat2-build ${1} ${1} 2>log.txt
    rm log.txt
fi
