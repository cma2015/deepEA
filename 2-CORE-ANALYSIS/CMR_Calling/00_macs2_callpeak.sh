echo "IP:${1}"
echo "input:${2}"
echo "SE or PE:${3}";
echo "genome size:${4}";
echo "Build model:${5}";
echo "p.value or q.value:${6}"
echo "measure:${7}"
echo "mfold_lower:${8}"
echo "mfold_upper:${9}"
echo "band width:${10}"
echo "extsize:${11}"
echo "shift:${12}"

if [ "${3}" == "BAM" ];then
bam="BAM"
else
bam='BAMPE'
fi

if [ "${6}" == "pvalue" ];then
measure="-p ${7} "
else
measure="-p ${7} "
fi
mkdir -p /home/galaxy/tools/2-CORE-ANALYSIS/CMR_Calling/out/
exp="macs2"
if [ "${5}" == "create_model" ];then
mfold="-m ${8} ${9} "
bw="--bw ${10} "
cmd="/home/miniconda2/bin/macs2 callpeak -t ${1} -c ${2} -g ${4} -f $bam $measure $mfold $bw --outdir /home/galaxy/tools/2-CORE-ANALYSIS/CMR_Calling/out/ -n $exp"
else
extsize="--extsize ${11} "
shift="--shift ${12} "
cmd="/home/miniconda2/bin/macs2 callpeak -t ${1} -c ${2} -g ${4} -f $bam $measure $extsize $shift --outdir /home/galaxy/tools/2-CORE-ANALYSIS/CMR_Calling/out/ --nomodel -n $exp"
fi

echo $cmd > /home/galaxy/tools/2-CORE-ANALYSIS/CMR_Calling/$$.sh && chmod u+x /home/galaxy/tools/2-CORE-ANALYSIS/CMR_Calling/$$.sh && bash /home/galaxy/tools/2-CORE-ANALYSIS/CMR_Calling/$$.sh && rm /home/galaxy/tools/2-CORE-ANALYSIS/CMR_Calling/$$.sh



