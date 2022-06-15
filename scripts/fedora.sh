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
dnf install python python2 pip
python2 -m ensurepip --no-default-pip

#Installing conda and pip dependencies
conda activate
conda install -y -c bioconda biopython matplotlib ipython jupyter pandas sympy nose seaborn psutil pysam sra-tools y kat busco=5.2.2 picard soapdenovo2 bwa bcftools trimmomatic samtools
pip3 install biopython matplotlib ipython jupyter pandas sympy nose seaborn psutil pysam

#Installing other dpendencies
dnf -y groupinstall "Development Tools" "Development Libraries"
dnf -y install samtools

mkdir dependencies
chmod 777 dependencies
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
#echo "Installing Redundans"
#dnf install wget curl git perl gcc g++ ldconfig
#pip2 install matplotlib numpy
#git clone --recursive https://github.com/lpryszcz/redundans.git
#cd redundans && bin/.compile.sh
#cd ..

echo "####################"
echo "Installing Redundans"
echo "####################"
conda env create -f $SELF/redundans_env.yml


echo "#################"
echo "Installing BUSCO"
echo "#################"
conda env create -f $SELF/busco_env.yml

#Installing nQuire
echo "Installing nQuire"
git clone --recursive https://github.com/clwgg/nQuire.git
chmod 777 nQuire
cd nQuire
make submodules
make
cd ..

#Installing Trimmomatic
echo "Installing Trimmomatic"
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$TRIMMOMATIC_VERSION.zip
unzip Trimmomatic-$TRIMMOMATIC_VERSION.zip

#Installing Picard-tools
echo "Installing Picard-Tools"
wget https://downloads.sourceforge.net/project/picard/picard-tools/$PICARD_VERSION/picard-tools-$PICARD_VERSION.zip
unzip picard-tools-$PICARD_VERSION.zip

cd ..


python3 "$SELF/../bin/create_config.py" --trimmomatic "$dep_folder/Trimmomatic-$TRIMMOMATIC_VERSION/" --karyon "$SELF/../" --redundans "$SELF/dependencies/redundans/" --SPAdes "$dep_folder/SPAdes-$SPAdes_VERSION-Linux" --nQuire "$dep_folder/nQuire/" --SOAPdenovo "$dep_folder/SOAPdenovo2-bin-LINUX-generic-r240" --redundans "$SELF/dependencies/redundans/" --picardtools "$dep_folder/picard-tools-$PICARD_VERSION" --output "$SELF/../configuration.txt"

