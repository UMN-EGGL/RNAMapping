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
    unzip \
    build-essential \
	apt-transport-https \
    openjdk-8-jre-headless \
    python3 \
    python3-dev \
    python3-pip \
    zlib1g-dev \
    libbz2-dev

RUN pip3 install ipython
RUN pip3 install snakemake
RUN pip3 install boto3
RUN wget -O adapterremoval-2.2.2.tar.gz https://github.com/MikkelSchubert/adapterremoval/archive/v2.2.2.tar.gz
RUN tar -xvzf adapterremoval-2.2.2.tar.gz
RUN cd adapterremoval-2.2.2; make
RUN wget https://github.com/alexdobin/STAR/archive/2.6.0a.tar.gz
RUN tar -xvzf 2.6.0a.tar.gz
RUN cd STAR-2.6.0a/source; make
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
RUN unzip fastqc_v0.11.7.zip
WORKDIR /RNAMapping
COPY . RNAMapping

# Build the Container with:
# $ docker build -t rnamap:latest .

# Run the Container passing through a port to the host
# $ docker run -p 4000:4000 -it rnamap

# Inside the container
# $ cd RNAMapping
