echo "reference genome:${1}"
echo "genome annotation:${2}"
echo "threads: ${3}"
echo "fastq: ${4}"
echo "method: ${5}"
echo "md: ${6}"
echo "cr: ${7}"
echo "mBQ: ${8}"
echo "mcov: ${9}"
echo "single-end or paired-end: ${10}"
echo "reversed fastq: ${11}"

toolDir="/home/galaxy/tools/DeepEA_software/meRanTK-1.2.0/"
outDir="/home/galaxy/tools/2-CORE-ANALYSIS/CMR_Calling/m5C/"
mkdir $outDir

if [ "${10}" == "SE" ];then
	if [ "${5}" == "meRanGs" ];then
		# generating a bisulfite index for meRanGs
		${toolDir}meRanGs mkbsidx -fa ${1} -GTF ${2} -t ${3} -GTFtagEPT Parent -GTFtagEPG gene -id ${outDir}meRanGsIDX

		# aligning RNA-BSseq reads to the test genome using the meRanGs index
		ln -sf ${4} ${4}.fastq
		${toolDir}meRanGs align -t ${3} -f ${4}.fastq -id ${outDir}meRanGsIDX -bg -o ${outDir} -S "meRanGs.sam" --star_genomeLoad NoSharedMemory -MM -un

		# calling m5Cs from meRanGs aligned RNA-BSseq reads
		${toolDir}meRanCall -p ${3} -s ${outDir}meRanGs_sorted.bam -f ${1} -gref -o ${outDir}meRanGs_meRanCall_m5Cs.txt -md ${6} -cr ${7} -mBQ ${8} -bed63 -sc 5 -mcov ${9}
	elif [ "${5}" == "meRanGh" ];then
		# generating a bisulfite index for meRanGh
		${toolDir}meRanGh mkbsidx  -t ${3} -fa ${1} -id ${outDir}meRanGhIDX

		# aligning RNA-BSseq reads to the test genome using the meRanGh index
		ln -sf ${4} ${4}.fastq
		${toolDir}meRanGh align -t ${3} -f ${4}.fastq  -id ${outDir}meRanGhIDX -GTF ${2} -bg -o ${outDir} -S "meRanGh.sam" -MM -un

		# calling m5Cs from meRanGh aligned RNA-BSseq reads
		${toolDir}meRanCall -p ${3} -s ${outDir}meRanGh_sorted.bam -f ${1} -gref -o ${outDir}meRanGh_meRanCall_m5Cs.txt -md ${6} -cr ${7} -mBQ ${8} -bed63 -sc 5 -mcov ${9}
	else
		# align RNA-BSseq reads to a set of reference transcripts
		# generating a refSeq transcript to genen name map file
		${toolDir}util/mkRefSeq2GeneMap.pl -f ${1} -m ${outDir}refSeq.t2g.map

		# generating a refSeq bisulfite index for meRanT
		${toolDir}meRanT mkbsidx -t ${3} -fa ${1} -id ${outDir}meRanTIDX

		# aligning RNA-BSseq reads to the test transcriptome using the meRanT index
		ln -sf ${4} ${4}.fastq
		${toolDir}meRanT align -t ${3} -f ${4}.fastq -i2g ${outDir}refSeq.t2g.map -o ${outDir} -S "meRanT.sam" -x ${outDir}refSeqRNA.500_C2T

		# calling m5Cs from meRanT aligned RNA-BSseq reads
		${toolDir}meRanCall -p ${3} -s ${outDir}meRanT_sorted.bam -f ${1} -tref -o ${outDir}meRanT_meRanCall_m5C.txt -md ${6} -cr ${7} -mBQ ${8} -sc 5 -mcov ${9}
	fi
else
	if [ "${5}" == "meRanGs" ];then
		# generating a bisulfite index for meRanGs
		${toolDir}meRanGs mkbsidx -fa ${1} -GTF ${2} -t ${3} -GTFtagEPT Parent -GTFtagEPG gene -id ${outDir}meRanGsIDX

		# aligning RNA-BSseq reads to the test genome using the meRanGs index
		ln -sf ${4} ${4}.fastq
		ln -sf ${11} ${11}.fastq
		${toolDir}meRanGs align -t ${3} -f ${4} -r ${11}.fastq -id ${outDir}meRanGsIDX -bg -o ${outDir} -S "meRanGs.sam" --star_genomeLoad NoSharedMemory -MM -un

		# calling m5Cs from meRanGs aligned RNA-BSseq reads
		${toolDir}meRanCall -p ${3} -s ${outDir}meRanGs_sorted.bam -f ${1} -gref -o ${outDir}meRanGs_meRanCall_m5Cs.txt -md ${6} -cr ${7} -mBQ ${8} -bed63 -sc 5 -mcov ${9}
	elif [ "${5}" == "meRanGh" ];then
		# generating a bisulfite index for meRanGh
		${toolDir}meRanGh mkbsidx  -t ${3} -fa ${1} -id ${outDir}meRanGhIDX

		# aligning RNA-BSseq reads to the test genome using the meRanGh index
		ln -sf ${4} ${4}.fastq
		ln -sf ${11} ${11}.fastq
		${toolDir}meRanGh align -t ${3} -f ${4}.fastq -r ${11}.fastq -id ${outDir}meRanGhIDX -GTF ${2} -bg -o ${outDir} -S "meRanGh.sam" -MM -un

		# calling m5Cs from meRanGh aligned RNA-BSseq reads
		${toolDir}meRanCall -p ${3} -s ${outDir}meRanGs_sorted.bam -f ${1} -gref -o ${outDir}meRanGh_meRanCall_m5Cs.txt -md ${6} -cr ${7} -mBQ ${8} -bed63 -sc 5 -mcov ${9}
	else
		# align RNA-BSseq reads to a set of reference transcripts
		# generating a refSeq transcript to genen name map file
		${toolDir}util/mkRefSeq2GeneMap.pl -f ${1}  -m ${outDir}refSeq.t2g.map

		# generating a refSeq bisulfite index for meRanT
		${toolDir}meRanT mkbsidx -t ${3} -fa ${1} -id ${outDir}meRanTIDX

		# aligning RNA-BSseq reads to the test transcriptome using the meRanT index
		ln -sf ${4} ${4}.fastq
		ln -sf ${11} ${11}.fastq
		${toolDir}meRanT align -t ${3} -f ${4}.fastq -r ${11}.fastq -i2g ${outDir}refSeq.t2g.map -o ${outDir} -S "meRanT.sam" -x ${outDir}refSeqRNA.500_C2T

		# calling m5Cs from meRanT aligned RNA-BSseq reads
		${toolDir}meRanCall -p ${3} -s ${outDir}meRanT_sorted.bam -f ${1} -tref -o ${outDir}meRanT_meRanCall_m5C.txt -md ${6} -cr ${7} -mBQ ${8} -sc 5 -mcov ${9}
	fi
fi

