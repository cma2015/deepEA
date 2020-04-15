FROM zhai2018/galaxy:19.05
MAINTAINER Jingjing Zhai, zhaijingjing603@gmail.com


ARG R_VERSION=3.6.1
ARG OS_IDENTIFIER=ubuntu-1804

# Install R
RUN apt-get install wget -y && \
    apt-get install curl -y && \
    apt-get update && \
    apt-get install libxml2 libxml2-dev -y && \
    apt-get update && \
    apt-get install libcurl4-openssl-dev libxml2-dev -y && \
    apt-get install libssl-dev -y && \
    apt install zlib1g -y && \
    apt install zlib1g-dev -y && \
    apt-get update -y && \
    apt-get install python-dev -y && \
    apt-get install libbz2-dev -y && \
    apt-get install -y liblzma-dev && \
    apt-get install libpng-dev -y && \
    apt-get install libjpeg-dev -y && \
    wget https://cdn.rstudio.com/r/${OS_IDENTIFIER}/pkgs/r-${R_VERSION}_1_amd64.deb && \
    apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -f -y ./r-${R_VERSION}_1_amd64.deb && \
    ln -s /opt/R/${R_VERSION}/bin/R /usr/bin/R && \
    ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/bin/Rscript && \
    ln -s /opt/R/${R_VERSION}/lib/R /usr/lib/R && \
    rm r-${R_VERSION}_1_amd64.deb && \
    rm -rf /var/lib/apt/lists/* && \
    curl -s -L https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh > ~/miniconda.sh \
    && /bin/bash ~/miniconda.sh -b -p /home/miniconda2/ \
    && rm ~/miniconda.sh && \
    echo "export PATH='/home/miniconda2/bin:$PATH'" >> ~/.bashrc && \
    /bin/bash -c "source ~/.bashrc" && \
    R -e "install.packages(c('BiocManager','remotes'), repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')" && \
    R -e "BiocManager::install(c('Rsamtools','rtracklayer','GenomicFeatures'))" && \
    wget http://bioconductor.statistik.tu-dortmund.de/packages/3.8/bioc/src/contrib/Guitar_1.20.1.tar.gz && \
    wget http://bioconductor.statistik.tu-dortmund.de/packages/3.0/bioc/src/contrib/exomePeak_1.6.0.tar.gz && \
    R -e "remotes::install_local(c('Guitar_1.20.1.tar.gz','exomePeak_1.6.0.tar.gz'))" && \
    rm Guitar_1.20.1.tar.gz exomePeak_1.6.0.tar.gz && \
    R -e "BiocManager::install(c('knitr', 'rmarkdown','highr', 'markdown', 'xfun', 'base64enc', 'tinytex'))" && \
    R -e "remotes::install_github(c('skgrange/threadr','skyhorsetomoon/Trumpet'))" && \
    git clone https://github.com/compgenomics/MeTPeak.git && \
    R CMD build MeTPeak && \
    R -e "install.packages('MeTPeak_1.0.0.tar.gz')" && \
    rm -r MeTPeak && rm MeTPeak_1.0.0.tar.gz && \
    R -e "BiocManager::install(c('DiffBind','plotly','RCAS','scales','gridExtra','diffloop','knitr','kableExtra','seqLogo','DiffLogo','seqinr','argparse','snowfall','data.table','BayesPeak','magrittr','BSgenome','GenomicScores', 'pipeR','caret','randomForest','DALEX','pheatmap','flexdashboard'))" && \
    mkdir -p /home/DeepEA/galaxy/tools/DeepEA_software && \
    cd /home/DeepEA/galaxy/tools/DeepEA_software && \
    wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz && \
    wget https://icbi.i-med.ac.at/software/meRanTK/downloads/1.2.0/meRanTK-1.2.0.zip && \
    tar -xzvf tophat-2.1.1.Linux_x86_64.tar.gz && rm tophat-2.1.1.Linux_x86_64.tar.gz && \
    unzip meRanTK-1.2.0.zip && rm -r meRanTK-1.2.0/testdata && \
    git clone https://github.com/ZW-xjtlu/m6ALogisticModel.git && \
    R -e "remotes::install_local('/home/DeepEA/galaxy/tools/DeepEA_software/m6ALogisticModel')" && \
    rm -r /home/DeepEA/galaxy/tools/DeepEA_software/m6ALogisticModel && \
    /home/miniconda2/bin/conda install -c bioconda/label/cf201901 sra-tools -y && \
    /home/miniconda2/bin/conda install -c bioconda fastp -y && \
    /home/miniconda2/bin/conda install -c bioconda fastqc -y && \
    /home/miniconda2/bin/conda install -c bioconda bowtie2 -y && \
    /home/miniconda2/bin/conda install -c bioconda star -y && \
    /home/miniconda2/bin/conda install -c bioconda hisat2 -y && \
    /home/miniconda2/bin/conda install -c bioconda bwa -y && \
    /home/miniconda2/bin/conda install -c bioconda meme -y && \
    /home/miniconda2/bin/conda install -c bioconda macs2 -y && \
    /home/miniconda2/bin/pip install biopython pysam numpy matplotlib scikit-learn scipy xgboost


# Copy xml
COPY 1-PRE-ANALYSIS /home/DeepEA/galaxy/tools/1-PRE-ANALYSIS
COPY 2-CORE-ANALYSIS /home/DeepEA/galaxy/tools/2-CORE-ANALYSIS
COPY 4-MULTI-OMICS_ANALYSIS /home/DeepEA/galaxy/tools/4-MULTI-OMICS_ANALYSIS
COPY galaxy.yml /home/DeepEA/galaxy/config/galaxy.yml
COPY tool_data_table_conf.xml /home/DeepEA/galaxy/config/tool_data_table_conf.xml
COPY tool_conf.xml /home/DeepEA/galaxy/config/tool_conf.xml
COPY welcome.html /home/DeepEA/galaxy/static/welcome.html
COPY assets /home/DeepEA/galaxy/static/assets
COPY demo_output /home/DeepEA/galaxy/static/demo_output
