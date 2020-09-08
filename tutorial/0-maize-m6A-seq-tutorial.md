## 0. Introduction
In this tutorial, we will show how to use deepEA to perform comprehensive m<sup>6</sup>A sequencing data analysis. The m<sup>6</sup>A -immunoprecipitated (IP) and input (non-IP) samples of *Zea mays* with two biological replicates were listed in the following table. More information regarding these maize m<sup>6</sup>A sequencing datasets is available in the reference (<a href="http://www.plantphysiol.org/content/182/1/332" target="_blank">Luo <I>et al</I>., 2020</a>).

|  Samples   | **Experiments** | Replicates  |
| :--------: | :-------------: | :---------: |
| SRR8383013 |       IP        | Replicate 1 |
| SRR8383014 |       IP        | Replicate 2 |
| SRR8383017 |      input      | Replicate 1 |
| SRR8383018 |      input      | Replicate 2 |
## 1. Download *Zea mays* reference genome sequences and annotation

Before analyzing m6A sequencing data, we firstly use the function **Obtain Genome Sequences and Annotation** in **Data Preparation** module to download *Zea mays* B73 reference genome sequences and GTF annotation, the following two screenshots shows details about how to execute this step:

#### Step 1: download *Zea mays* reference genome sequences in FASTA format

![1](../assets/img/maize-m6A-seq-tutorial/1.jpeg)

#### Step 2: download *Zea mays* genome annotation file in GTF format

![2](../assets/img/maize-m6A-seq-tutorial/2.jpeg)

After that, reference genome sequences (named as **zea_mays_release-47_Genome.txt**) and annotation (named as **zea_mays_release-47_GTF.txt**) will be listed in your History Panel.

## 2. Download raw m<sup>6</sup>A sequencing reads

After finishing downloading *Zea mays* B73 reference genome sequences and annotation, we start to download raw m<sup>6</sup>A sequencing reads from NCBI SRA's database, this process can be finished by the function **Obtain Epitranscriptome Sequencing Reads** in **Data Preparation** module. Please see the following screenshots for details:

#### Step 1: download IP sample for biological replicate 1 (SRR8383013)

![3](../assets/img/maize-m6A-seq-tutorial/3.jpeg)

#### Step 2: download IP sample for biological replicate 2 (SRR8383014)
![4](../assets/img/maize-m6A-seq-tutorial/4.jpeg)

#### Step 3: download input sample for biological replicate 1 (SRR8383017)
![5](../assets/img/maize-m6A-seq-tutorial/5.jpeg)

#### Step 4: download input sample for biological replicate 2 (SRR8383018)
![6](../assets/img/maize-m6A-seq-tutorial/6.jpeg)

After finishing the above four steps, four raw m<sup>6</sup>A sequencing reads named as **SRR8383013.sra, SRR8383013.sra, SRR8383013.sra, SRR8383013.sra**, respectively will be listed in you History Panel

## 3. Convert raw m<sup>6</sup>A sequencing reads in SRA to FASTQ format

Here, let's convert raw m<sup>6</sup>A sequencing reads from SRA to FASTQ format by the function **Sequencing Data Preprocessing** in **Data Preparation** module, please see the following four screenshots for details:

#### Step 1: process SRR8383013.sra

![7](../assets/img/maize-m6A-seq-tutorial/7.jpeg)

#### Step 2: process SRR8383014.sra
![8](../assets/img/maize-m6A-seq-tutorial/8.jpeg)

#### Step 3: process SRR8383017.sra
![9](../assets/img/maize-m6A-seq-tutorial/9.jpeg)

#### Step 4: process SRR8383018.sra
![10](../assets/img/maize-m6A-seq-tutorial/10.jpeg)

For each SRA accession, two FASTQ format files (forward reads and reverse reads) will be generated as this experiment is paired-end sequencing.

## 4. Trim raw m<sup>6</sup>A sequencing reads

Before align m<sup>6</sup>A sequencing reads to genome, trim low-quality reads is necessary in all NGS (Next Generation Sequencing) analyses as which may cause incorrect mapping of reads to a reference genome, and even result in incorrect identification of RNA modifications. deepEA provided the function **Assess Reads Quality** in **Quality Control** module to filter raw reads to clean reads, please see the following screenshots for details:

#### Step 1: trim SRR8383013

![11](../assets/img/maize-m6A-seq-tutorial/11.jpeg)

#### Step 2: trim SRR8383014
![12](../assets/img/maize-m6A-seq-tutorial/12.jpeg)

#### Step 3: trim SRR8383017
![13](../assets/img/maize-m6A-seq-tutorial/12.jpeg)

#### Step 4: trim SRR8383018
![14](../assets/img/maize-m6A-seq-tutorial/14.jpeg)

For each SRA accession, three files will be output, please see the following screenshot for detail:

![27](../assets/img/maize-m6A-seq-tutorial/27.jpeg)

## 5. Align clean m<sup>6</sup>A sequencing reads to *Zea mays* B73 reference genome with HISAT2

Here, we start to align clean reads to reference genome with HISAT2 provided in module **Identification of RNA Modifications**, see the following screenshots for details:

#### Step 1: align SRR8383013

![15](../assets/img/maize-m6A-seq-tutorial/15.jpeg)

#### Step 2: align SRR8383014
![16](../assets/img/maize-m6A-seq-tutorial/16.jpeg)

#### Step 3: align SRR8383017
![17](../assets/img/maize-m6A-seq-tutorial/17.jpeg)

#### Step 4: align SRR8383018
![18](../assets/img/maize-m6A-seq-tutorial/18.jpeg)

Then for each SRA accession, reads-genome alignments in BAM format and alignment summary will be generated, see the following screenshot:
![28](../assets/img/maize-m6A-seq-tutorial/28.jpeg)

## 6. Call m<sup>6</sup>A enriched peaks with macs2

After finishing aligning reads to genome, let's start to call m<sup>6</sup>A enriched peaks with macs2, the following screenshots show details about parameter settings:

#### Step 1: call m<sup>6</sup>A peaks for biological replicate 1

![19](../assets/img/maize-m6A-seq-tutorial/19.jpeg)

This step will generate ~16,400 peaks for biological replicate 1

#### Step 2: call m<sup>6</sup>A peaks for biological replicate 2

![20](../assets/img/maize-m6A-seq-tutorial/20.jpeg)

This step will generate ~16,100 peaks for biological replicate 2

#### Step 3: obtain consistent peaks between two biological replciates

After finishing peak calling for two replicates, deepEA also provided a function **Merge two biological replicates** to obtain consistent peaks between two biological replicates, see the following screenshots for details:

![21](../assets/img/maize-m6A-seq-tutorial/21.jpeg)

Then the consistent peaks named as **intersect.bed** (about ~14,000 consistent peaks) will be shown in you History Panel.

## 7. Perform functional annotation for m<sup>6</sup>A
#### m<sup>6</sup>A distribution

The distribution of m<sup>6</sup>A  in the genome and transcriptome can be visualized by the function **RNA Modification Distribution** in **Functional Annotation** module, see the folloing screenshot for detail:

![29](../assets/img/maize-m6A-seq-tutorial/29.jpeg)

Then an interactively HTML will be generated, please see <a href="https://cma2015.github.io/deepEA/tutorial/demo/RNA_modifications_distribution.html" target="_blank">here</a> to preview this results

#### Link m<sup>6</sup>A modifications with gene
![22](../assets/img/maize-m6A-seq-tutorial/22.jpeg)

This output for this function are shown in the following the screenshot
![30](../assets/img/maize-m6A-seq-tutorial/30.jpeg)

#### *De-novo* motif discovery
This following screenshot shows how to use homer to perform *de-novo* motif discovery
![32](../assets/img/maize-m6A-seq-tutorial/32.jpeg)

Then an HTML document will be generated, please click <a href="https://cma2015.github.io/deepEA/tutorial/demo/homer_motifs.html" target="_blank">here</a> to preview this results

#### GO functional enrichment analysis
deepEA provided the function **Functional Enrichment Analysis** to perform GO enrichment analysis, see the following screenshot for detail:

![31](../assets/img/maize-m6A-seq-tutorial/31.jpeg)

Then a figure in PDF format and a TAB seperated matrix will be generated, please click  <a href="https://cma2015.github.io/deepEA/tutorial/demo/GO_enrichment.pdf" target="_blank">here</a>   for details.

## 7. Multi-omics integrative analysis

To run this module, you have to download test data provided by deepEA, and then upload the data in directory `test_data/Multi-omics Integrative Analysis/` to deepEA server. If you are not sure how to upload local data into deepEA server, please see <a href="https://www.youtube.com/watch?reload=9&v=vDd9yQHiYYQ" target="_blank">here</a> for details. Then you can run the function **Integrative Analysis of Three Omics Data Sets** in **Multi-omics Integrative Analysis** module as the following screenshot shows:

![33](../assets/img/maize-m6A-seq-tutorial/33.jpeg)

Then an interactively HTML document will be output, click <a href="https://cma2015.github.io/deepEA/tutorial/demo/multi-omics.html" target="_blank">here</a>  to preview this results.

## 8. Build an m<sup>6</sup>A predictor based on machine learning

#### Step 1: generate positive and negative samples for m<sup>6</sup>A predictor construction
![23](../assets/img/maize-m6A-seq-tutorial/23.jpeg)
#### Step 2: encoding positve samples
![24](../assets/img/maize-m6A-seq-tutorial/24.jpeg)
#### Step 3: encoding negative samples
![25](../assets/img/maize-m6A-seq-tutorial/25.jpeg)
#### Step 4: m<sup>6</sup>A predictor construction and evaluation
![26](../assets/img/maize-m6A-seq-tutorial/26.jpeg)