#!/bin/python
import sys, os, re, subprocess, math
import argparse
import psutil
import pysam
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
	parser.add_argument('-o', '--output_name', required=True, help="Output prefix")
	parser.add_argument('-v', '--vcf', required=True, help="VCF file used as input")
	parser.add_argument('-p', '--pileup', required=True, help="Mpileup file used as input")
	parser.add_argument('-b', '--bam', required=True, help="Bam file used as input")
	parser.add_argument('-l', '--library', required=True, help="Illumina library used for the KAT plot")
	parser.add_argument('--configuration', default=False, help="Configuration file. By default will use ./configuration.txt as the configuration file.")
	parser.add_argument('-w', '--wsize', default=1000, help="Window size for plotting")
	parser.add_argument('-x', '--max_scaf2plot', default=20, help="Number of scaffolds to analyze")
	parser.add_argument('-s', '--scafminsize', default=False, help="Will ignore scaffolds with length below the given threshold")
	parser.add_argument('-S', '--scafmaxsize', default=False, help="Will ignore scaffolds with length above the given threshold")
	parser.add_argument('-i', '--job_id', default=False, help='Identifier of the intermediate files generated by the different programs. If false, the program will assign a name consisting of a string of 6 random alphanumeric characters.')

args = parser.parse_args()


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

def id_generator(size=6, chars=string.ascii_uppercase + string.digits): 
	return ''.join(random.choice(chars) for _ in range(size))
	
config_path = args.configuration
if not args.configuration:
	selfpath = os.path.dirname(os.path.realpath(sys.argv[0]))
	config_path = selfpath[:selfpath.rfind('/')]
	config_path = selfpath[:selfpath.rfind('/')]+"/configuration.txt"
config_dict = parse_config(config_path)

	
counter = int(args.max_scaf2plot)
window_size=int(args.wsize)
step=window_size/2

true_output = os.path.abspath(args.output_directory)
cwd = os.path.abspath(os.getcwd())
os.chdir(true_output)

os.system("bgzip -c "+ args.vcf + " > " + args.vcf + ".gz")
os.system("tabix -p vcf "+ args.vcf+".gz")
#vcf_file = pysam.VariantFile(args.vcf+".gz", 'r')
bam_file = pysam.AlignmentFile(args.bam, 'rb')
home = config_dict["karyon"][0]
job_ID = args.job_id if args.job_id else id_generator()
name = args.output_name if args.output_name else job_ID
kitchen = home + "kitchen/"+job_ID

lendict = {}
fastainput = SeqIO.index(args.fasta, "fasta")
for i in fastainput:
	lendict[i] = len(fastainput.get_raw(i).decode())

from karyonplots import katplot, allplots
allplots(window_size, 
				args.vcf, 
				args.fasta, 
				args.bam, 
				args.pileup, 
				args.library, 
				config_dict['nQuire'][0], 
				config_dict["KAT"][0], 
				kitchen, 
				true_output, 
				counter, 
				job_ID, name,
	 			args.scafminsize,
	 			args.scafmaxsize
				)

os.chdir(cwd)	



	
	
	
