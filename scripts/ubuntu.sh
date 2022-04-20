#!/bin/sh

SELF="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
mkdir ../tmp

GATK_VERSION=4.1.9.0
SPAdes_VERSION=3.9.0
HTSLIB_VERSION=1.9
SAMTOOLS_VERSION=1.9
BCFTOOLS_VERSION=1.9
BWA_VERSION=0.7.15
TRIMMOMATIC_VERSION=0.36
PICARD_VERSION=1.78

#Install python and pip
apt-get install python python2 pip
python2 -m ensurepip --no-default-pip

#Installing conda and pip dependencies
conda install -c bioconda -y biopython matplotlib ipython jupyter pandas sympy nose seaborn psutil pysam sra-tools gatk4 kat busco=5.2.2 picard soapdenovo2 bwa bcftools trimmomatic
pip3 install biopython matplotlib ipython jupyter pandas sympy nose seaborn psutil pysam

#Installing other dpendencies
apt-get update
apt-get install -y libbz2-dev
apt-get install -y bzip2
apt-get install -y zlib1g-dev
apt-get install -y libncurses5-dev 
apt-get install -y libncursesw5-dev
apt-get install -y liblzma-dev
apt-get -y install samtools

mkdir dependencies
cd dependencies
dep_folder=`pwd`

#Installins SPADes
SPAdes_VERSION=3.9.0
wget https://github.com/ablab/spades/releases/download/v$SPAdes_VERSION/SPAdes-$SPAdes_VERSION-Linux.tar.gz
tar -xvf SPAdes-$SPAdes_VERSION-Linux.tar.gz

#Installing SOAPdeNovo
echo "Installing SOAPdeNovo"
wget https://downloads.sourceforge.net/project/soapdenovo2/SOAPdenovo2/bin/r240/SOAPdenovo2-bin-LINUX-generic-r240.tgz
tar -xvf SOAPdenovo2-bin-LINUX-generic-r240.tgz

#Installing Redundans
apt-get install wget curl git perl gcc g++ ldconfig
pip2 install matplotlib numpy
git clone --recursive https://github.com/lpryszcz/redundans.git
cd redundans && bin/.compile.sh
cd ..

#Installing nQuire
git clone --recursive https://github.com/clwgg/nQuire.git
chmod 777 nQuire
cd nQuire
make submodules
make
cd ..

#Installing Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$TRIMMOMATIC_VERSION.zip
unzip Trimmomatic-$TRIMMOMATIC_VERSION.zip

cd ..


python3 "$SELF/../bin/create_config.py" --trimmomatic "$dep_folder/Trimmomatic-$TRIMMOMATIC_VERSION/" --karyon "$SELF/../" --redundans "$SELF/dependencies/redundans/" --SPAdes "$dep_folder/SPAdes-$SPAdes_VERSION-Linux" --nQuire "$dep_folder/nQuire/" --SOAPdenovo "$dep_folder/SOAPdenovo2-bin-LINUX-generic-r240" --redundans "$SELF/dependencies/redundans/" --output "$SELF/../configuration.txt"

