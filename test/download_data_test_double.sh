#!/bin/bash
###
# Karyon Test .
# version 0.1a
###
echo "Starting sratoolkit download..."
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
./sratoolkit.2.10.6-ubuntu64/bin/vdb-config --interactive
echo "Starting data download..."
./sratoolkit.2.10.6-ubuntu64/bin/fastq-dump -split-files DRR040667
