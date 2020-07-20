#!/bin/bash
###
# Karyon dependencies installer for UNIX.
# version 0.1a
###

echo -n "Installing dependencies...\n"
apt-get update
apt-get install -y python
apt-get install -y git
apt-get install -y wget
apt-get install -y unzip
apt-get install -y nano
apt-get install -y gcc
apt-get install -y g++
apt-get install -y make
apt-get install -y unzip
apt-get install -y python3-pip
apt-get install -y wget libxml-libxml-perl
apt-get install -y default-jre
apt-get install -y tabix
apt-get install -y autoconf automake libtool
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz -O /tmp/sratoolkit.tar.gz \
	&& tar zxvf /tmp/sratoolkit.tar.gz -C /opt/ && rm /tmp/sratoolkit.tar.gz
PATH="/opt/sratoolkit.2.10.0-ubuntu64/bin/:${PATH}"
