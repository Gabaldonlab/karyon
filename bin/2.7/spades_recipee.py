#!/bin/python
import sys, numpy, os.path, re
import argparse
from Bio import SeqIO

def call_SPAdes(library_file, path, output, name, commands, no_diploid, memory_limit, nodes):
	libstring = ' '
	backstring = ''
	chunk_col = 5
	for i in open(library_file):
		chunk = i.split()		
		if chunk[chunk_col] == "1":
			libstring = libstring + " -1 " + os.path.abspath(chunk[0]) + " -2 " + os.path.abspath(chunk[6]) + " "
		elif chunk[chunk_col] == "2": 
			continue
		elif chunk[chunk_col] == "--12":
			libstring = libstring + " --12 " + os.path.abspath(chunk[0])
		else:
			backstring = backstring + os.path.abspath(chunk[chunk_col]) + " " + os.path.abspath(chunk[0]) + " "
	
	libstring = libstring + backstring
	outputfile = open(output+name+"_karyon.job", 'w')
	if no_diploid == True:
		outputfile.write("python2 " + os.path.join(path, "bin", "spades.py") + libstring + " -t " + str(nodes) + " -m " +  str(memory_limit) + " " +commands + "-o " + os.path.join(output, "spades"))
	else:
		outputfile.write("python2 " + os.path.join(path, "bin", "dipspades.py") + libstring + " -t " + str(nodes) + " -m " +  str(memory_limit) + " " +commands + "-o " + output)
	outputfile.write("\n")
	outputfile.close()
