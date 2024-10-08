# Karyon
![Latest Version](https://img.shields.io/github/v/tag/gabaldonlab/karyon?label=Latest%20Version)
![Docker Pulls](https://img.shields.io/docker/pulls/cgenomics/karyon)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://hub.docker.com/repository/docker/cgenomics/karyon)

This repository contains the Karyon pipeline.

## Introduction

Karyon is a pipeline for the assembly and analysis of highly heterozygous genomes. It uses redundans (Pryszcz & Gabaldón, 2016) to reduce heterozygosity during the assembly process, and then maps the original libraries against the reduced assembly to analyze the distribution of heterozygous regions. With this information, it generates a series of plots that can aid researchers to generate informed hypotheses with regard of the architecture of their genomes.

Scripts contained in this repository:
1) **karyon.py** -> complete pipeline, including genome assembly, assembly reduction, SNP calling and plot generation
2) **prepare_libraries.py** -> karyon dependency. It uses Trimmomatic () to trim input libraries before genome assembly.
3) **spades_recipee.py** -> Karyon dependency. It generates a file that launches dipSPAdes () with the input.
4) **varcall_recipee.py** -> Karyon dependency. It generates a file that launches all steps in the SNP calling pipeline.
5) **karyonplots.py** -> Karyon dependency. It generates all the plots as part of the Karyon pipeline.
6) **allplots.py** -> Standalone version of karyonplots.py. It allows the user to input karyon results to generate the plots again.
7) **nQuire_plot.py** -> It allows the user to run the local ploidy plot alone.
8) **Dockerfile** -> Docker file involved in building the image to run Karyon
9) **install.sh** -> Bash script required to install the remaining dependencies in the dockerfile
10) **redundans_env.yml & busco_env.yml** -> Conda environments in YAML format required to install some of the trickiest dependencies.

## How to install it

* **Quick start**
You can install it using the standard installation or through Docker.

1. **Standard installation**

Follow this steps


```Shell
# First clone the Karyon repository
git clonehttps://github.com/Gabaldonlab/karyon.git
# Change to karyon/scripts directory
cd karyon/scripts/
# Then, run the installation script.
bash installation.sh
```
2. **Docker installation**

In order to run this container you'll need docker installed. Need to get [started](https://docs.docker.com/get-started/)?

- [Windows](https://docs.docker.com/desktop/windows/install/)
- [OS X](https://docs.docker.com/desktop/mac/install/)
- [Linux](https://docs.docker.com/desktop/linux/install/)
  
* **Use the Dockerfile build**

```Shell
# From the karyon git directory
docker build --no-cache -t cgenomics/karyon:1.2 .
# Start the container and indicate a volume and a container name
docker run -dit --name=karyon -v $(pwd):/root/src/karyon/shared --rm cgenomics/karyon:1.2
# Install all the necessary dependencies inside the running container. First run interactively the container
docker exec -it karyon bash
# Changing dir to the karyon volume in the container where the Dockerfile is located
cd /root/src/karyon/shared/
# Run the dependency installation script
bash scripts/docker_install.sh
```

* **Pull the docker image from Docker Hub**

Pull `gabaldonlab/karyon` from the Docker repository:
```Shell
# First pull the image
docker pull cgenomics/karyon:1.2
#  Start the container and indicate a volume and a container name
docker run -dit --name=karyon -v $(pwd):/root/src/karyon/shared --rm cgenomics/karyon:1.2
# Install all the necessary dependencies inside the running container. First run interactively the container
docker exec -it karyon bash
# Changing dir to the karyon volume in the container where the Dockerfile is located
cd /root/src/karyon/shared/
# Run the dependency installation script
bash scripts/docker_install.sh
```
## Test dataset

The test dataset is composed by two sequencing libraries from NCBI SRA corresponding to Lichtheimia ramosa B5399, one of the strains analyzed in the main publication. 
```Shell
# Execute interactively the docker container
docker exec -it karyon bash
# Configure SRA tools within the docker container
~/src/karyon/shared/dependencies/sratoolkit.3.0.0-ubuntu64/bin/vdb-config --interactive
# Download the SRA libraries at the desired location
cd /root/src/karyon/shared/
~/src/karyon/shared/dependencies/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump --split-files SRR974799 SRR974800

```

## Manual

Please, check the [manual](https://github.com/Gabaldonlab/karyon/blob/master/Karyon_manual.pdf) for a comprehensive use of Karyon

## Authors 
* **Miguel Ángel Naranjo Ortiz** - *Pipeline work* - [MANaranjo](https://github.com/MANaranjo)
* **Manuel Molina Marín** - *Docker work* - [manumolina](https://github.com/manumolina)
* **Diego Fuentes Palacios** - *Docker work & testing* [dfupa](https://github.com/dfupa)
* **Toni Gabaldón** - *Intellectual design & validation* - [tgabaldon](https://github.com/tgabaldon)

## License 
This project is licensed under the GNU General Public License - see the [LICENSE.md](LICENSE.md) file for details.
