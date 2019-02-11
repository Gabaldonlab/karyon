#!/bin/bash
###
# Karyon installer for UNIX.
# bash <(curl -Ls https://github.com/Gabaldonlab/karyon/karyon_installer)
# version 0.1a
###

branch="master" 
if [ ! -z $1 ]; then branch=$1; fi

echo "#####################################################################################"
echo "#                                                                                   #"
echo "#                                   Karyon installer                                #"
echo "#                                                                                   #"
echo "#        version 0.1a                                  Miguel Angel Naranjo         #"
echo "#####################################################################################"
echo ""
echo "Karyon and its dependencies will be installed in:" `pwd`/karyon
echo "Installation will take 5-10 minutes. "
# echo "To track the installation status execute in the new terminal:"
# echo "  tail -f `pwd`/karyon/$log"
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

# real software to install
# awk wget tar git gcc g++ make cd ln date ldconfig unzip python samtools

error=""
# check if all basic programs exists
for cmd in echo awk git wget unzip tar nano gcc g++ make cd ln date ldconfig unzip pip; do
    if ! exists $cmd; then
        case $cmd in        
            "pip")
                echo "Install pip first (ie. 'sudo apt-get install python-pip')!"  
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
    echo -e "\nAborted due to missing dependencies (see above)!"
    return 1;
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
git clone https://github.com/clwgg/nQuire.git

echo "Installing anaconda_ete"
echo "To-Do"

echo "Installing Python packages"
pip install --upgrade pip
pip install numpy
pip install biopython
pip install psutil
pip install pysam
python -m pip install --user matplotlib ipython jupyter pandas sympy nose seaborn

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


PATH="$(pwd)/samtools-1.9/:${PATH}"
PATH="$(pwd)/bcftools-1.9/:${PATH}"

echo "Installing BWA"
wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2
tar -vxjf bwa-0.7.15.tar.bz2
cd bwa-0.7.15
make 
cd ..

PATH="$(pwd)/bwa-0.7.15/:${PATH}"

apt-get clean
set -x; rm -rf /var/lib/apt/lists/*

# git clone -b $branch --recursive https://github.com/Gabaldonlab/karyon.git >> /dev/null 2>&1 
# cd karyon

# compile dependencies
# sh bin/.compile.sh `pwd`/$log
# retcode=$?; 
# if [ $retcode -gt 0 ]; then
#     echo "  [$retcode] ERROR!"
#     tail -n 20 $log
#     return $retcode
# fi


echo 'alias karyon="python /root/src/karyon/scripts/karyon.py"' >> ~/.bashrc


echo `date` "Installation finished!"
echo ""
echo "To try Karyon, execute:"
# To-Do
# echo "cd karyon; ./karyon.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1"
# echo ""
# echo "To uninstall execute:"
# echo "rm -rI `pwd`"
# echo ""
# To-Do
echo "#####################################################################################"
echo "# Karyon depends on several programs (http://bit.ly/redundans_dependencies)      #"
echo "# Acknowledge their authors, get familiar with licensing and register if necessary. #"
echo "#####################################################################################"
echo ""