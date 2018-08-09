FROM ubuntu:18.04
MAINTAINER Sam<beeso018@umn.edu>
LABEL Description "https://github.com/UMN-EGGL/RNAMapping"

# Install the necessary packages ontop of base ubuntu installation 
RUN apt-get -y update && apt-get install -y \
    curl \
	lsb-release \ 
    wget \
    git \
    gcc \
    build-essential \
	apt-transport-https \
    python3 \
    python3-dev \
    python3-pip \
    zlib1g-dev \
    libbz2-dev

RUN cd /root

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -f
ENV PATH=/root/miniconda3/bin:${PATH}
RUN conda update -n base conda

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN conda create -y -n default python=3 fastqc star adapterremoval
RUN /bin/bash -c "source activate default"

RUN pip3 install ipython
RUN pip3 install snakemake
RUN pip3 install boto3

COPY . root/RNAMapping

# Build the Container with:
# $ docker build -t rnamap:latest .

# Run the Container passing through a port to the host
# $ docker run -p 4000:4000 -it rnamap

# Inside the container
# $ cd /root/RNAMapping
# $ source activate default
