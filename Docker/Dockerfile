FROM ubuntu:18.04
MAINTAINER mbrown@broadinstitute.org
ENV DEBIAN_FRONTEND=noninteractive

#~~~~~~~~~~~~~~
# Install Preliminaries
#~~~~~~~~~~~~~~
RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends \
    automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
    g++ \
    gcc \
    git \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libjpeg-dev \
    liblzma-dev \
    libpng-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    make \
    pbzip2 \
    pigz \
    python \
    python-pip \
    python-setuptools \
    python3-pip \
    r-base \
    unzip \
    wget \
    x11-common \
    zlib1g \
    zlib1g-dev \
    zlibc

# make python3 be the default python
# RUN ln -sf /usr/bin/python3 /usr/bin/python


#~~~~~~~~~~~
# Python 
#~~~~~~~~~~~
RUN pip install --upgrade pip && pip install \
	argparse==1.4.0 \
    cutadapt

RUN pip3 install --upgrade pip && pip3 install \
    argparse==1.4.0 \
    cutadapt

#~~~~~~~~~~~
# Environmental Variables 
#~~~~~~~~~~~
ENV SRC /usr/local/src
ENV BIN /usr/local/bin



#~~~~~~~~~~~
## Samtools
#~~~~~~~~~~~
WORKDIR $SRC
ENV SAMTOOLS_VERSION=1.12
RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    cd $SRC && \
    wget $SAMTOOLS_URL && \
    tar xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
    cd ../ && ./configure --without-curses && make && make install && \
    rm $SRC/samtools-${SAMTOOLS_VERSION}.tar.bz2 && rm -r $SRC/samtools-${SAMTOOLS_VERSION}


#~~~~~~~~~~~
## Bedtools
#~~~~~~~~~~~
RUN wget -q https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
   tar xf bedtools-2.30.0.tar.gz && \
   cd bedtools2 && \
   make && \
   cp bin/* $BIN/ && \
   rm $SRC/bedtools-2.30.0.tar.gz && rm -r $SRC/bedtools2

#~~~~~~~~~~~~~~
# BWA
#~~~~~~~~~~~~~~
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    cp -r * $BIN/ && \
    cd .. && \
    rm -r bwa

#~~~~~~~~~~~~~~
# FastQC
#~~~~~~~~~~~~~~
WORKDIR $SRC
ENV FASTQC_VERSION=0.11.9
RUN FASTQC_URL="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" && \
    cd $SRC && \
    wget $FASTQC_URL && \
    unzip fastqc_v${FASTQC_VERSION}.zip && \
    cp -r FastQC $BIN && \
    rm -r FastQC




#~~~~~~~~~~~~~~
# TopHat2
#~~~~~~~~~~~~~~

WORKDIR $SRC
ENV TopHat2_VERSION=2.1.1
# RUN wget https://ccb.jhu.edu/software/tophat/downloads/tophat-${TopHat2_VERSION}.OSX_x86_64.tar.gz && \
RUN wget http://ccb.jhu.edu/software/tophat/downloads/tophat-${TopHat2_VERSION}.Linux_x86_64.tar.gz && \
    tar zxvf tophat-${TopHat2_VERSION}.Linux_x86_64.tar.gz && \
    cd $BIN && \
    ln -s $SRC/tophat-${TopHat2_VERSION}.Linux_x86_64/tophat2 . 



#~~~~~~~~~~~~~~
# BowTie2
#~~~~~~~~~~~~~~
WORKDIR $SRC
ENV Bowtie2_VERSION=2.4.4

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-source.zip && \
    unzip bowtie2-2.4.4-source.zip && \
    cd bowtie2-2.4.4 && \
    make && \
    cp -r * $BIN


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## VirTect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

WORKDIR $SRC
ARG CACHEBUST=2
RUN git clone https://github.com/brownmp/VirTect.git 


# RUN wget -q https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh -O $SRC/CTAT-VirusIntegrationFinder/WDL/monitor_script.sh
# RUN chmod a+rx $SRC/CTAT-VirusIntegrationFinder/WDL/monitor_script.sh


RUN apt-get -qq -y remove automake build-essential cmake curl git gcc g++ wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log


