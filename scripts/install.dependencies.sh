#!/bin/bash
###
# Karyon dependencies installer for UNIX.
# version 0.1a
###

echo -n "Installing dependencies...\n"
apt-get update
apt-get install -y git
apt-get install -y wget
apt-get install -y unzip
apt-get install -y nano
apt-get install -y gcc
apt-get install -y g++
apt-get install -y make
apt-get install -y unzip
apt-get install -y python-pip