echo "genome: ${1}"
if [ ! -d "${1}.star" ]; then
    mkdir -p ${1}".star"
    STAR --runMode genomeGenerate --genomeDir  ${1}".star" --genomeSAindexNbases 10 --genomeFastaFiles ${1} --runThreadN 4;
fi
