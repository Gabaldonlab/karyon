FROM python:3.9.1

# metadata
LABEL base.image="python:3.9.1"
LABEL version="1"
LABEL software="Karyon"
LABEL software.version="1.1"
LABEL description=""
LABEL website="https://github.com/Gabaldonlab/karyon"
LABEL license="GNU General Public License"
LABEL maintainer="Manu Molina (BSC), Diego Fuentes (BSC)"

#Set up miniconda
ENV PATH="/root/miniconda3/bin:$PATH"
ARG PATH="/root/miniconda3/bin:$PATH"
RUN apt-get update
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

SHELL ["/bin/bash", "-c"]
RUN conda --version

#Install dependencies for the main environment in python3
RUN apt-get update
RUN apt-get install -y software-properties-common
RUN apt-get install -y build-essential
RUN apt-get -y install git curl nano zip sshfs autoconf automake libtool apt-utils gcc g++ make libxml-libxml-perl default-jre tabix unzip
RUN apt-get install -y python3-pip python3-lxml libapache2-mod-wsgi-py3
RUN pip3 install --upgrade ipython jupyter numpy biopython python-dev-tools pandas sympy nose seaborn psutil pysam
