#!/bin/bash
###
# Karyon installer for UNIX.
# version 0.1a
###

branch="master" 
if [ ! -z $1 ]; then branch=$1; fi

echo "##################################################################################################"
echo "#                                                                                                #"
echo "#                                         Karyon installer                                       #"
echo "#                                                                                                #"
echo "#        version 0.1a                                             Miguel Angel Naranjo           #"
echo "##################################################################################################"
echo ""
echo "Karyon and its dependencies will be installed in:" `pwd`/karyon
echo "Installation will take 5-10 minutes. "
echo ""

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

'''
Dependencies versions
'''
GATK_VERSION=4.0.12.0
SPAdes_VERSION=3.9.0
HTSLIB_VERSION=1.9
SAMTOOLS_VERSION=1.9
BCFTOOLS_VERSION=1.9
BWA_VERSION=0.7.15
TRIMMOMATIC_VERSION=0.36
PICARD_VERSION=1.78


echo " Creating dependencies folder..."
mkdir dependencies
rm -r ../kitchen
mkdir ../kitchen
cd dependencies

echo " Installing basic software..."
apt-get install -y software-properties-common
apt-get install -y build-essential

echo " Installing Java..."
add-apt-repository -y ppa:webupd8team/java
apt-get update
echo debconf shared/accepted-oracle-license-v1-1 select true | debconf-set-selections
echo debconf shared/accepted-oracle-license-v1-1 seen true | debconf-set-selections
apt-get install -y --force-yes oracle-java8-installer

echo "Installing GATK..."
wget https://github.com/broadinstitute/gatk/releases/download/$GATK_VERSION/gatk-$GATK_VERSION.zip
unzip gatk-$GATK_VERSION.zip

echo "Installing SOAPdeNovo"
wget https://downloads.sourceforge.net/project/soapdenovo2/SOAPdenovo2/bin/r240/SOAPdenovo2-bin-LINUX-generic-r240.tgz
tar -xvf SOAPdenovo2-bin-LINUX-generic-r240.tgz

echo "Installing SPAdes"
wget https://github.com/ablab/spades/releases/download/v$SPAdes_VERSION/SPAdes-$SPAdes_VERSION-Linux.tar.gz
tar -xvf SPAdes-$SPAdes_VERSION-Linux.tar.gz

echo "Installing Trimmomatic"
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$TRIMMOMATIC_VERSION.zip
unzip Trimmomatic-$TRIMMOMATIC_VERSION.zip

echo "Installing Picard-tools"
wget https://downloads.sourceforge.net/project/picard/picard-tools/$PICARD_VERSION/picard-tools-$PICARD_VERSION.zip
unzip picard-tools-$PICARD_VERSION.zip

echo "Installing nQuire"
git clone --recursive https://github.com/clwgg/nQuire.git
sudo chmod 777 nQuire
cd nQuire
make submodules
make
cd ..

echo "Installing anaconda_ete"
conda install kat
#wget https://repo.continuum.io/miniconda/Miniconda3-3.7.0-Linux-x86_64.sh -O ~/miniconda.sh
#bash ~/miniconda.sh -b -p $HOME/miniconda

echo "Installing Python packages"
# pip install --upgrade pip
python3 -m pip install --upgrade pip
pip3 install numpy
pip3 install biopython
pip3 install psutil
pip3 install pysam
python3 -m pip install --user matplotlib ipython jupyter pandas sympy nose seaborn

echo "Installing Samtools, Bcftools and Htslib..."
apt-get update
apt-get install -y libbz2-dev
apt-get install -y bzip2
apt-get install -y zlib1g-dev
apt-get install -y libncurses5-dev 
apt-get install -y libncursesw5-dev
apt-get install -y liblzma-dev

wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
tar -vxjf htslib-1.9.tar.bz2
cd htslib-1.9
make
cd ..

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make
cd ..

wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9
make
cd ..



echo "Installing BWA"
wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2
tar -vxjf bwa-0.7.15.tar.bz2
cd bwa-0.7.15
make 
cd ..

echo "Installing Redundans"
git clone --recursive https://github.com/lpryszcz/redundans.git
cd redundans && bin/.compile.sh
cd ..


rm ./*.bz2
rm ./*.zip
rm ./*.zip.1
rm ./*.tar.gz
rm ./*.tgz

dep_folder=/home/dependencies
PATH="$HOME/miniconda/bin:${PATH}"
PATH="$dep_folder/samtools-1.9/:${PATH}"
PATH="$dep_folder/bcftools-1.9/:${PATH}"
PATH="$dep_folder/bwa-0.7.15/:${PATH}"
echo 'alias karyon="python $(pwd)/bin/2.7/karyon.py"' >> ~/.bashrc

cd ..
chmod -R 777 dependencies
cd dependencies

apt-get clean
set -x; rm -rf /var/lib/apt/lists/*

echo `date` "Installation finished!"
echo "##################################################################################################"
echo "# Karyon depends on several programs (https://github.com/Gabaldonlab/karyon#prerequisities)      #"
echo "# Acknowledge their authors, get familiar with licensing and register if necessary.              #"
echo "##################################################################################################"
echo ""
