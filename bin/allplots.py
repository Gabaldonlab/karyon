#!/bin/python
import sys, os, re, subprocess, math
import argparse
import psutil
from pysam import pysam
from Bio import SeqIO
import numpy as np
import numpy.random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import seaborn as sns
import pandas as pd
import scipy.stats
from scipy.stats import gaussian_kde
from scipy import stats
from decimal import Decimal
import string, random


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--fasta', required=True, help="fasta file used as input")
	parser.add_argument('-d', '--output_directory', default="./", help='Directory where all the output files will be generated.')
	parser.add_argument('-o', '--output_name', required=False, help="Output prefix")
	parser.add_argument('-v', '--vcf', required=True, help="VCF file used as input")
	parser.add_argument('-p', '--pileup', required=True, help="Mpileup file used as input")
	parser.add_argument('-b', '--bam', required=True, help="Bam file used as input")
	parser.add_argument('-l', '--library', required=True, nargs='+',  help="Illumina libraries used for the KAT plot")
	parser.add_argument('--configuration', default=False, help="Configuration file. By default will use ./configuration.txt as the configuration file.")
	parser.add_argument('-w', '--window_size', default=1000, help="Window size for plotting")
	parser.add_argument('-x', '--max_scaf2plot', default=20, help="Number of scaffolds to analyze")
	parser.add_argument('-s', '--scafminsize', default=False, help="Will ignore scaffolds with length below the given threshold")
	parser.add_argument('-S', '--scafmaxsize', default=False, help="Will ignore scaffolds with length above the given threshold")
	parser.add_argument('-i', '--job_id', default=False, help='Identifier of the intermediate files generated by the different programs. If false, the program will assign a name consisting of a string of 6 random alphanumeric characters.')

args = parser.parse_args()
true_output = os.path.abspath(args.output_directory)
if true_output[-1] != "/":
	true_output=true_output+"/"

def parse_config(config):
	config_dict = {}
	prev = 0
	for line in open(config):
		if line[0] == "#": continue
		elif line[0] == "+":
			prev = line[1:-1]
			config_dict[prev] = ["","",""]
		elif line[0] == "@":
			if config_dict[prev][0] != "": continue
			config_dict[prev][0] = line[1:-1]
		elif line[0] == ">":
			config_dict[prev][1] = config_dict[prev][1] + line[1:-1] + " "
		elif line[0] == "?":
			if config_dict[prev][2] != "": continue
			config_dict[prev][2] = line[1:-1] + " "
	return config_dict

def id_generator(fastainput):
	return(fastainput[fastainput.rfind("/")+1:fastainput.rfind(".")])	 
	
config_path = args.configuration
if not args.configuration:
	selfpath = os.path.dirname(os.path.realpath(sys.argv[0]))
	config_path = selfpath[:selfpath.rfind('/')]
	config_path = selfpath[:selfpath.rfind('/')]+"/configuration.txt"
config_dict = parse_config(config_path)

	
counter = int(args.max_scaf2plot)
window_size=int(args.window_size)
step=window_size/2

cwd = os.path.abspath(os.getcwd())
os.chdir(true_output)

os.system("bgzip -c "+ args.vcf + " > " + args.vcf + ".gz")
os.system("tabix -p vcf "+ args.vcf+".gz")
bam_file = pysam.AlignmentFile(args.bam, 'rb')
home = config_dict["karyon"][0]
job_ID = args.job_id if args.job_id else id_generator(args.fasta)
name = args.output_name if args.output_name else job_ID
name = name[name.rfind("/")+1:]
kitchen = home + "tmp/"+job_ID

lendict = {}
fastainput = SeqIO.index(args.fasta, "fasta")
for i in fastainput:
	lendict[i] = len(fastainput[i].seq)

from karyonplots import katplot, allplots
from report import report, ploidy_veredict
df = allplots(window_size, 
				args.vcf, 
				args.fasta, 
				args.bam, 
				args.pileup, 
				args.library[0], 
				config_dict['nQuire'][0], 
				config_dict["KAT"][0], 
				kitchen, 
				true_output, 
				counter, 
				job_ID, name,
	 			args.scafminsize,
	 			args.scafmaxsize, False)

df2 = ploidy_veredict(df, true_output, name, window_size)
report(true_output, name, df2, True, False, window_size, False, False)
df2.to_csv(true_output+"/Report/"+name+".csv", index=False)
os.chdir(cwd)	



	
	
	
