#!/bin/bash
###
# Karyon dependencies installer for Docker.
# version 1.1
###

branch="master" 
if [ ! -z $1 ]; then branch=$1; fi

echo "##################################################################################################"
echo "#                                                                                                #"
echo "#                                         Karyon installer                                       #"
echo "#                                                                                                #"
echo "#        version 1.1                      Miguel Angel Naranjo & Diego Fuentes Palacios          #"
echo "##################################################################################################"
echo ""
echo "Karyon and its dependencies will be installed in:" `pwd`/karyon
echo "Installation will take 5-10 minutes. "
echo ""

SELF="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
mkdir ../tmp
echo "$SELF"
# sleep
echo "I'll proceed with installation in 2 seconds... Press Ctrl-C to cancel."
sleep 2s

#echo "Changing the python interpreter"
#rm -f /usr/bin/python3 && ln -s /usr/bin/python3.6 /usr/bin/python3
#echo "Done"

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
                echo "Install pip first (ie. 'sudo apt-get install python-pip3')!"  
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
        echo "Installing"
        bash $SELF/install.dependencies.sh
    fi
fi

echo " Everything looks good :) Let's proceed..."
sleep 2s

pwd=`pwd`
init_path="$(pwd)";


#Dependencies versions

GATK_VERSION=4.1.9.0
SPAdes_VERSION=3.9.0
HTSLIB_VERSION=1.9
SAMTOOLS_VERSION=1.9
BCFTOOLS_VERSION=1.9
BWA_VERSION=0.7.15
TRIMMOMATIC_VERSION=0.36
PICARD_VERSION=1.78


echo " Creating dependencies folder..."
mkdir $SELF/../dependencies
cd $SELF/../dependencies

echo "####################"
echo "Installing SRA tools"
echo "####################"

wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.0/sratoolkit.3.0.0-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
export PATH=$PATH:$PWD/sratoolkit.3.0.0-ubuntu64/bin >> ~/.bashrc

echo "####################"
echo "Installing Redundans"
echo "####################"
conda env create -f $SELF/redundans_env.yml


echo "#################"
echo "Installing BUSCO"
echo "#################"
conda env create -f $SELF/busco_env.yml

echo "####################"
echo "Installing Platanus"
echo "####################"
mkdir platanus-1.2.4
mkdir platanus-1.2.4/bin
cd platanus-1.2.4/
cd bin/
wget -O- http://platanus.bio.titech.ac.jp/?ddownload=145 > platanus && chmod +x platanus
cd $SELF/../dependencies
export PATH=$PATH:$PWD/platanus-1.2.4/bin/ >> ~/.bashrc

#echo "###############"
#echo "Installing KAT"
#echo "###############"
#git clone https://github.com/TGAC/KAT.git
#cd KAT
#./build_boost.sh
#./autogen.sh
#./configure
#make install

echo "###################"
echo " Installing Java..."
echo "###################"
add-apt-repository -y ppa:webupd8team/java
apt-get update
echo debconf shared/accepted-oracle-license-v1-1 select true | debconf-set-selections
echo debconf shared/accepted-oracle-license-v1-1 seen true | debconf-set-selections
apt-get install -y --force-yes oracle-java8-installer

echo "###################"
echo "Installing GATK..."
echo "###################"
wget https://github.com/broadinstitute/gatk/releases/download/$GATK_VERSION/gatk-$GATK_VERSION.zip
unzip gatk-$GATK_VERSION.zip

echo "#####################"
echo "Installing SOAPdeNovo"
echo "#####################"
wget https://downloads.sourceforge.net/project/soapdenovo2/SOAPdenovo2/bin/r240/SOAPdenovo2-bin-LINUX-generic-r240.tgz
tar -xvf SOAPdenovo2-bin-LINUX-generic-r240.tgz

echo "##################"
echo "Installing SPAdes"
echo "##################"
wget https://github.com/ablab/spades/releases/download/v$SPAdes_VERSION/SPAdes-$SPAdes_VERSION-Linux.tar.gz
tar -xvf SPAdes-$SPAdes_VERSION-Linux.tar.gz
export PATH=$PATH:$PWD/SPAdes-$SPAdes_VERSION-Linux/bin >> ~/.bashrc
echo "######################"
echo "Installing Trimmomatic"
echo "######################"
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-$TRIMMOMATIC_VERSION.zip
unzip Trimmomatic-$TRIMMOMATIC_VERSION.zip

echo "#######################"
echo "Installing Picard-tools"
echo "#######################"
wget https://downloads.sourceforge.net/project/picard/picard-tools/$PICARD_VERSION/picard-tools-$PICARD_VERSION.zip
unzip picard-tools-$PICARD_VERSION.zip

echo "##################"
echo "Installing nQuire"
echo "##################"
git clone --recursive https://github.com/clwgg/nQuire.git
chmod 777 nQuire
cd nQuire
make submodules
make
cd ..


echo "###########################################"
echo "Installing Samtools, Bcftools and Htslib..."
echo "###########################################"
apt-get update
apt-get install -y libz-dev
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

echo "###############"
echo "Installing BWA"
echo "###############"
wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2
tar -vxjf bwa-0.7.15.tar.bz2
cd bwa-0.7.15
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
#PATH="$CONDASITE:${PATH}"
export PATH=$PWD/samtools-1.9/:${PATH} >> ~/.bashrc
export PATH=$PWD/bcftools-1.9/:${PATH} >> ~/.bashrc
export PATH=$PWD/bwa-0.7.15/:${PATH} >> ~/.bashrc

apt-get clean
set -x; rm -rf /var/lib/apt/lists/*

python3 "$SELF/../bin/create_config.py" --karyon "$SELF/../" --BWA "$dep_folder/bwa-0.7.15/" --samtools "$dep_folder/samtools-1.9/" --bcftools "$dep_folder/bcftools-1.9/" --picardtools "$dep_folder/picard-tools-$PICARD_VERSION" --GATK "$dep_folder/gatk-$GATK_VERSION" --SPAdes "$dep_folder/SPAdes-$SPAdes_VERSION-Linux" --nQuire "$dep_folder/nQuire/" --SOAPdenovo "$dep_folder/SOAPdenovo2-bin-LINUX-generic-r240" --trimmomatic "$dep_folder/Trimmomatic-$TRIMMOMATIC_VERSION/" --Platanus "$dep_folder/platanus-1.2.4/" --output "$SELF/../configuration.txt"

echo `date` "Installation finished!"
echo "##################################################################################################"
echo "# Karyon depends on several programs (https://github.com/Gabaldonlab/karyon#prerequisities)      #"
echo "# Acknowledge their authors, get familiar with licensing and register if necessary.              #"
echo "##################################################################################################"
echo ""
