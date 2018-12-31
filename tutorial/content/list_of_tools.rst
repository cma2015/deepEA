The tools
=========

.. note:: With the release of deepTools 2.0, we renamed a couple of tools:

    * **heatmapper** to :doc:`tools/plotHeatmap`
    * **profiler** to :doc:`tools/plotProfile`
    * **bamCorrelate** to :doc:`tools/multiBamSummary`
    * **bigwigCorrelate** to :doc:`tools/multiBigwigSummary`
    * **bamFingerprint** to :doc:`tools/plotFingerprint`.

 For more, see :doc:`changelog`.

.. contents:: 
    :local:

+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
| tool                                | type             | input files                         | main output file(s)                        | application                                                                       |
+=====================================+==================+=====================================+============================================+===================================================================================+
|:doc:`tools/multiBamSummary`         | data integration | 2 or more BAM                       | interval-based table of values             | perform cross-sample analyses of read counts --> plotCorrelation, plotPCA         |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/multiBigwigSummary`      | data integration | 2 or more bigWig                    | interval-based table of values             |  perform cross-sample analyses of genome-wide scores --> plotCorrelation, plotPCA |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotCorrelation`         | visualization    | bam/multiBigwigSummary output       | clustered heatmap                          | visualize the Pearson/Spearman correlation                                        |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotPCA`                 | visualization    | bam/multiBigwigSummary output       | 2 PCA plots                                | visualize the principal component analysis                                        |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotFingerprint`         | QC               | 2 BAM                               | 1 diagnostic plot                          | assess enrichment strength of a ChIP sample                                       |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/computeGCBias`           | QC               | 1 BAM                               | 2 diagnostic plots                         | calculate the exp. and obs. GC distribution of reads                              |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/correctGCBias`           | QC               | 1 BAM, output from computeGCbias    | 1 GC-corrected BAM                         | obtain a BAM file with reads distributed according to the genome’s GC content     |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/bamCoverage`             | normalization    | BAM                                 | bedGraph or bigWig                         | obtain the normalized read coverage of a single BAM file                          |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/bamCompare`              | normalization    | 2 BAM                               | bedGraph or bigWig                         | normalize 2 files to each other (e.g. log2ratio, difference)                      |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/computeMatrix`           | data integration | 1 or more bigWig, 1 or more BED     | zipped file for plotHeatmap or plotProfile | compute the values needed for heatmaps and summary plots                          |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/estimateReadFiltering`   | information      | 1 or more BAM files                 | table of values                            | estimate the number of reads filtered from a BAM file or files                    |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/alignmentSieve`          | QC               | 1 BAM file                          | 1 filtered BAM or BEDPE file               | filters a BAM file based on one or more criteria                                  |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotHeatmap`             | visualization    | computeMatrix output                | heatmap of read coverages                  | visualize the read coverages for genomic regions                                  |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotProfile`             | visualization    | computeMatrix output                | summary plot (“meta-profile”)              | visualize the average read coverages over a group of genomic regions              |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotCoverage`            | visualization    | 1 or more BAM                       | 2 diagnostic plots                         | visualize the average read coverages over sampled genomic  positions              |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/bamPEFragmentSize`       | information      | 1  BAM                              | text with paired-end fragment length       | obtain the average fragment length from paired ends                               |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/plotEnrichment`          | visualization    | 1 or more BAM and 1 or more BED/GTF | A diagnostic plot                          | plots the fraction of alignments overlapping the given features                   |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+
|:doc:`tools/computeMatrixOperations` | miscellaneous    | 1 or more BAM and 1 or more BED/GTF | A diagnostic plot                          | plots the fraction of alignments overlapping the given features                   |
+-------------------------------------+------------------+-------------------------------------+--------------------------------------------+-----------------------------------------------------------------------------------+

