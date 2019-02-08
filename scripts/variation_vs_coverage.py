#!/bin/python

import sys
import matplotlib.pyplot as plt
import numpy as np
import numpy.random
import argparse
from Bio import SeqIO
from scipy.stats import gaussian_kde

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', required=True, help="VCF file used as input")
    parser.add_argument('-p', '--pileup', required=True, help="Mpileup file used as input")
    parser.add_argument('-w', '--wsize', required=True, help="Window size for plotting")
    parser.add_argument('-f', '--fasta', default=False, help="Fasta file used as reference genome. Only used if scafsize is set")
    parser.add_argument('-s', '--scafsize', default=False, help="Will ignore scaffolds with length below the given threshold")
    parser.add_argument('-S', '--scafmaxsize', default=False, help="Will ignore scaffolds with length above the given threshold")

args = parser.parse_args()

vcf_file = open(args.vcf)
pileup_file = open(args.pileup)
window_size = int(args.wsize)

if args.scafsize != False:
	scaf_size = int(args.scafsize)
	fastainput = SeqIO.index(args.fasta, "fasta")
if args.scafmaxsize != False:
	scaf_max_size = int(args.scafmaxsize)
	fastainput = SeqIO.index(args.fasta, "fasta")

quality_filter = 500

lendict = {}
if args.fasta != False:
	for i in fastainput:
		lendict[i] = len(fastainput.get_raw(i).decode())


def extract_vcf_data (vcf_file):
	snp_count = 0
	curr_pos = 0
	curr_scaffold = ''
	snp_dict = {}
	for line in vcf_file:
		if line[0] == "#": continue
		chunk = line.split("\t")
		if args.scafsize != False:
			if (lendict[chunk[0]]) <= scaf_size:continue
		if args.scafmaxsize != False:
			if (lendict[chunk[0]]) >= scaf_max_size: continue
		if float(chunk[5]) < quality_filter: continue
		if chunk[0] != curr_scaffold:
			curr_scaffold = chunk[0]
			curr_pos = 0
			snp_dict[curr_scaffold+"_"+str(curr_pos)] = snp_count
			snp_count = 1
		elif int(chunk[1]) >= curr_pos + window_size:
			snp_dict[curr_scaffold+"_"+str(curr_pos)] = snp_count
			if (int(chunk[1])-curr_pos)/window_size > 1:
				for i in range(1, (int(chunk[1])-curr_pos)/window_size):
					snp_dict[curr_scaffold+"_"+str(curr_pos+(i*window_size))] = 0
			curr_pos = curr_pos + (window_size*((int(chunk[1])-curr_pos)/window_size))
			snp_count = 1
		else:
			snp_count = snp_count + 1
	return snp_dict

def extract_pileup_data (pileup_file):
	coverage = []
	curr_pos = 0
	curr_scaffold = ''
	cov_dict = {}
	for line in pileup_file:
		if line[0] == "#": continue
		chunk = line.split()
		if args.scafsize != False:
			if (lendict[chunk[0]]) <= scaf_size:continue
		if args.scafmaxsize != False:
			if (lendict[chunk[0]]) >= scaf_max_size: continue
		if chunk[0] != curr_scaffold:
			curr_scaffold = chunk[0]
			curr_pos = 0
			cov_dict[curr_scaffold+"_"+str(curr_pos)] = numpy.mean(coverage)
			coverage = []
		elif int(chunk[1]) >= curr_pos + window_size:
			cov_dict[curr_scaffold+"_"+str(curr_pos)] = numpy.mean(coverage)
			if (int(chunk[1])-curr_pos)/window_size > 1:
				for i in range(1, (int(chunk[1])-curr_pos)/window_size):
					cov_dict[curr_scaffold+"_"+str(curr_pos+(i*window_size))] = 0
			curr_pos = curr_pos + (window_size*((int(chunk[1])-curr_pos)/window_size))
			coverage = [int(chunk[3])]
		else: 
			coverage.append(int(chunk[3]))
	return cov_dict
	

snp_density = extract_vcf_data(vcf_file)
mean_cov = extract_pileup_data(pileup_file)

x, y = [], []
for element in snp_density:
	if element in mean_cov and mean_cov[element] < 500 and mean_cov[element] > 10 and snp_density[element] < 500:
		x.append(snp_density[element])
		y.append(mean_cov[element])

xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)
plt.figure(figsize=(3,2.5))
fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=10, edgecolor='')


plt.savefig('var_v_cov.png')





