import pysam
from Bio import SeqIO
import argparse
import numpy as np


def stopRatio(pulldownBAM, inputBAM, cDNA, out, stopNumber, stopRate, difference):
     pulldownSAM = pysam.AlignmentFile(pulldownBAM, "rb")
     inputSAM = pysam.AlignmentFile(inputBAM, "rb")
     transSeq = SeqIO.parse(cDNA, "fasta")
     outRatio = open(out, 'w+')
     for fasta in transSeq:
          name, sequence = fasta.id, str(fasta.seq)
          idx = [pos for pos, char in enumerate(sequence) if char == "A" or char == "a"]
          pulldownPosVec = [tt.pos for tt in pulldownSAM.fetch(name, 0, len(sequence))]
          inputPosVec = [tt.pos for tt in inputSAM.fetch(name, 0, len(sequence))]
          if len(pulldownPosVec) == 0 & len(inputPosVec) == 0:
               continue
          for j in idx:
               # the number of reads reading through position i in pulldown and input samples
               pull_readthroughtCount = pulldownSAM.count(name, j, j + 1)
               input_readthroughtCount = inputSAM.count(name, j, j + 1)

               # the number of reads with the mapping position starting at base i+1 in pulldown and input samples
               pull_stopCount = len([pos for pos, char in enumerate(pulldownPosVec) if char == j + 2])
               input_stopCount = len([pos for pos, char in enumerate(inputPosVec) if char == j + 2])
               if (pull_stopCount + pull_readthroughtCount) == 0:
                    continue
               else:
                   pull_stopRate = float(pull_stopCount)/(float(pull_stopCount) + float(pull_readthroughtCount))

                   if (input_stopCount + input_readthroughtCount) == 0:
                       input_stopRate = 0
                   else:
                        input_stopRate = float(input_stopCount) / (float(input_stopCount) + float(input_readthroughtCount))
                   outRatio.write("\t".join([name, str(j + 1), str(pull_readthroughtCount), str(pull_stopCount), str(pull_stopRate),
                                              str(input_readthroughtCount), str(input_stopCount), str(input_stopRate)]))
                   outRatio.write('\n')
     outRatio.close()
     inputSAM.close()
     pulldownSAM.close()
     results = np.loadtxt(fname=out, delimiter='\t', dtype='string')
     results = results[np.where(results[:, 3].astype('int') >= stopNumber)]
     results = results[np.where(results[:, 7].astype('float') <= stopRate)]
     results = results[np.where((results[:, 4].astype('float') - results[:, 7].astype('float')) >= difference)]
     np.savetxt(fname=out, delimiter='\t', X=results, header="\t".join(["transcripts", 'position', 'pulldown.readthrought.count',
                               'pulldown.stop.count', 'pulldown.stoprate','input.readthrought.count',
                               'input.stop.count', 'input.stoprate']), fmt="%s")
     return results


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", dest = "pulldown", type = str, default = None,
                        help = "The directory of pulldown bam.")
    parser.add_argument("-i", dest ="input", type = str, default = None,
                        help = "The directory of input bam.")
    parser.add_argument("-c", dest="cDNA", type=str, default=None,
                        help="The directory of cDNA sequence.")
    parser.add_argument("-out", dest = "outDir", help = "The directory of output")
    parser.add_argument("-s", dest="stopNumber", type=int, default = 5,
                        help = "The least number of reads in the N3-CMC(+) sample.")
    parser.add_argument("-r", dest="stopRate", type=float, default=0.1,
                        help="The least number of reads in the N3-CMC(-) sample.")
    parser.add_argument("-d", dest="difference", type=float, default=0.3,
                        help="The difference of stop rate between N3-CMC(+) and N3-CMC(-) sample.")
    # parser.add_argument("-f", dest="fold", type=int, default=2,
    #                     help="the difference of stop rate of position i must be at least twofold "
    #                          "greater than that of the two flanking bases.")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    res = stopRatio(pulldownBAM=args.pulldown, inputBAM=args.input, cDNA=args.cDNA,
                    out=args.outDir + "/stop_ratio.txt", stopNumber=args.stopNumber,
                    stopRate=args.stopRate, difference=args.difference)

