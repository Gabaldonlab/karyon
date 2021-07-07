#!/bin/bash
###
# Karyon installer for UNIX.
# version 0.1a
###

branch="master" 
if [ ! -z $1 ]; then branch=$1; fi

echo "##################################################################################################"
echo "#                                                                                                #"
echo "#                                   Karyon light installer                                       #"
echo "#                                                                                                #"
echo "#        version 0.1a                                             Miguel Angel Naranjo           #"
echo "##################################################################################################"
echo ""
echo "This will only install dependencies for allplots.py"
echo "For a full installation please run install.sh"
echo "Karyon and its dependencies will be installed in:" `pwd`/karyon
echo "Installation will take 5-10 minutes. "
echo ""

SELF="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
mkdir ../tmp

# sleep
echo "I'll proceed with installation in 2 seconds... Press Ctrl-C to cancel."
sleep 2s

echo ""
echo `date` "Checking dependencies..."

exists()
{
  command -v "$1" >/dev/null 2>&1
}

error=""
# check if all basic programs exists
for cmd in echo awk git wget unzip tar nano gcc g++ make cd ln date ldconfig unzip pip; do
    if ! exists $cmd; then
        case $cmd in        
            "pip")
                echo "Install pip first (ie. 'sudo apt-get install python-pip2')!"  
                ;;
            *)
                echo "Install $cmd first (ie. 'sudo apt-get install $cmd')!"
                ;;
        esac
        error=1
    fi
done

# check if all libs present #BWA
for lib in libz; do
    if [ -z "$(ldconfig -p | grep $lib.so)" ] ; then
        echo " Missing library $lib !"
        error=1
    fi
done

# check headers #BWA
for lib in zlib.h; do
    if [ ! -s /usr/include/$lib ] && [ ! -s /usr/lib/$lib ]; then
        echo " Missing headers $lib !"
        error=1
    fi
done

# skip if error
if [ ! -z $error ]; then
    echo -n "We've finded missing dependencies.\n"
    echo -n "Would you like to install (y or n)?\n"
    read response
    if [ $response != "y" ]; then        
        echo "\nAborted due to missing dependencies (see above)!"
        return 1;
    else
        sh install.dependencies.sh
    fi
fi

# check python version

PyVer=`python --version 2>&1 | cut -f2 -d" " | cut -f-2 -d"."`
if [ $PyVer != "2.7" ] && [ $PyVer != "2.6" ]; then 
    echo ""
    echo "[ERROR] Install Python 2.7!"
    echo "If you have Python 2.7 already installed, you can either "
    echo "make an alias before installation and use of Redundans ('alias python=python2.7' should do)"
    echo "or use Python virtual environment (https://virtualenv.pypa.io)."
    return 1
fi
echo " Everything looks good :) Let's proceed..."
sleep 2s

pwd=`pwd`
init_path="$(pwd)";

'''
Dependencies versions
'''
GATK_VERSION=4.1.9.0
SPAdes_VERSION=3.9.0
HTSLIB_VERSION=1.9
SAMTOOLS_VERSION=1.9
BCFTOOLS_VERSION=1.9
BWA_VERSION=0.7.15
TRIMMOMATIC_VERSION=0.36
PICARD_VERSION=1.78


echo " Creating dependencies folder..."
mkdir dependencies
cd $SELF/dependencies

echo " Installing basic software..."
apt-get install -y software-properties-common
apt-get install -y build-essential

CONDACHECK=`conda list | wc -l`
if [ $CONDACHECK < 2 ]; then
    echo "Installing Bioconda"
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    bash ~/miniconda.sh -b -p ~/miniconda 
    rm ~/miniconda.sh
    export PATH="$PATH:/root/miniconda/bin"
    echo 'alias conda="/root/miniconda/bin/conda"' >> ~/.bashrc
    source ~/.bashrc
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    echo "Bioconda OK"
    sleep 2s
    echo "Installing Bioconda"
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    bash ~/miniconda.sh -b -p ~/miniconda 
    rm ~/miniconda.sh
    export PATH="$PATH:/root/miniconda/bin"
    echo 'alias conda="/root/miniconda/bin/conda"' >> ~/.bashrc
    source ~/.bashrc
    echo "Bioconda OK"
    sleep 2s
            fi
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

echo "Installing Python packages"
conda install -y biopython 
conda install -y matplotlib 
conda install -y ipython
conda install -y jupyter
conda install -y pandas
conda install -y sympy
conda install -y nose
conda install -y seaborn
conda install -y psutil
conda install -y pysam
conda install -c bioconda -y sra-tools
pip3 install  biopython matplotlib ipython jupyter pandas sympy nose seaborn psutil pysam

echo "Installing KAT"
conda install -y kat
CONDAPATH=`which conda`
CONDASITE=$(echo "$CONDAPATH" | sed "s/\/conda//")
echo $CONDASITE
echo "alias kat=$CONDASITE/kat" >> ~/.bashrc
source ~/.bashrc

echo " Installing Java..."
add-apt-repository -y ppa:webupd8team/java
apt-get update
echo debconf shared/accepted-oracle-license-v1-1 select true | debconf-set-selections
echo debconf shared/accepted-oracle-license-v1-1 seen true | debconf-set-selections
apt-get install -y --force-yes oracle-java8-installer

echo "Installing nQuire"
git clone --recursive https://github.com/clwgg/nQuire.git
chmod 777 nQuire
cd nQuire
make submodules
make
cd ..

rm ./*.bz2
rm ./*.zip
rm ./*.zip.1
rm ./*.tar.gz
rm ./*.tgz

cd ..
chmod -R 777 dependencies
cd dependencies

dep_folder=`pwd`
PATH="$dep_folder/samtools-1.9/:${PATH}"
PATH="$dep_folder/bcftools-1.9/:${PATH}"
PATH="$dep_folder/bwa-0.7.15/:${PATH}"

apt-get clean
set -x; rm -rf /var/lib/apt/lists/*

python3 "$SELF/../bin/create_config.py" --karyon "$SELF/../" --samtools "$dep_folder/samtools-1.9/" --nQuire "$dep_folder/nQuire/" --output "$SELF/../configuration.txt"

echo `date` "Installation finished!"
echo "##################################################################################################"
echo "# Karyon depends on several programs (https://github.com/Gabaldonlab/karyon#prerequisities)      #"
echo "# Acknowledge their authors, get familiar with licensing and register if necessary.              #"
echo "##################################################################################################"
echo ""
