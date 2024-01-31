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
import pandas as pd
import scipy.stats
from scipy.stats import gaussian_kde
from scipy import stats
from decimal import Decimal
import string, random

parser = argparse.ArgumentParser(description="")
parser.add_argument("-c", "--configuration", default=False, help="Configuration file. By default will use ./configuration.txt as the configuration file.")
parser.add_argument("-Q", "--nQuire", default=False, help="Location of nQuire. If ommited, it will look for the program from the configuration file")
parser.add_argument("-d", "--output_directory", required=True, help="Directory to send all the output files")
parser.add_argument('-t', '--temporary', default="./tmp", help="Directory to dump temporary files.")
parser.add_argument("-v", "--vcf_file", required=True, help="VCF file, required as input.")
parser.add_argument("-f", "--fasta_file", required=True, help="Reference fasta file, required as input.")
parser.add_argument("-b", "--bam_file", required=True, help="Bam file, required as input.")
parser.add_argument("-w", "--window_size", default=1000, help="Size of the sliding windows to analyze. Default is 1000bp.")
parser.add_argument('-s', '--scafminsize', default=False, help="Will ignore scaffolds with length below the given threshold")
parser.add_argument('-S', '--scafmaxsize', default=False, help="Will ignore scaffolds with length above the given threshold")
parser.add_argument('-x', '--max_scaf2plot', default=20, help="Maximum number of scaffolds to plot for scaffold-specific plots. Default is 20.")

args = parser.parse_args()

window_size, step = int(args.window_size), int(args.window_size)/2

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


if args.nQuire != False:
	nQuire = os.abspath(args.nQuire)+"/nQuire"
else:
	if args.configuration != False:
		config_path = os.abspath(args.configuration)
	else:
		selfpath = os.path.dirname(os.path.realpath(sys.argv[0]))
		config_path = selfpath[:selfpath.rfind('/')]+"/configuration.txt"
	config_dict = parse_config(config_path)
	nQuire = config_dict["nQuire"][0]
	
true_output = os.path.abspath(args.output_directory)
if true_output[-1] != "/":
	true_output=true_output+"/"
if args.temporary[-1] != "/":
	args.temporary=args.temporary+"/"

config_dict = parse_config(config_path)
home = config_dict["karyon"][0]
if home[-1] != "/": home = home + "/"

fastainput = SeqIO.index(args.fasta_file, "fasta")
lendict = {}
for i in fastainput:
	lendict[i] = len(fastainput.get_raw(i).decode())

from karyonplots import window_walker
window_walker(window_size, step, args.vcf_file, args.fasta_file, args.bam_file, nQuire, args.temporary, true_output, int(args.max_scaf2plot), lendict, int(args.scafminsize), int(args.scafmaxsize), True)
