#!/bin/python
import sys, os, re, subprocess, math
import argparse
import psutil
import pysam
from Bio import SeqIO
import numpy as np
import numpy.random
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns
import pandas as pd
import scipy.stats
from decimal import Decimal
from scipy import stats



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-o', '--output', required=True, help="Output prefix")
	parser.add_argument('-v', '--vcf', required=True, help="VCF file used as input")
	parser.add_argument('-p', '--pileup', required=True, help="Mpileup file used as input")
	parser.add_argument('-w', '--wsize', default=1000, help="Window size for plotting")
	parser.add_argument('-s', '--scafsize', default=False, help="Will ignore scaffolds with length below the given threshold")
	parser.add_argument('-S', '--scafmaxsize', default=False, help="Will ignore scaffolds with length above the given threshold")
	
args = parser.parse_args()
window_size = int(args.wsize)

def extract_pileup_data (pileup_file):
	coverage = []
	curr_pos = 0
	curr_scaffold = ''
	cov_dict = {}
	scaf_size=int(args.scafsize)
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
	pileup_file.seek(0)
	return cov_dict

def extract_vcf_data (vcf_file, quality_filter, scaf):
	snp_count = 0
	curr_pos = 0
	curr_scaffold = ''
	snp_dict = {}
	scaf_size = int(args.scafsize)
	for line in vcf_file:
		if line[0] == "#": continue
		if scaf != '' and line.split()[0] != scaf: continue
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
	vcf_file.seek(0)
	return snp_dict

def var_v_cov(vcf, pileup):
	vcf_file = open(vcf)
	pileup_file = open(pileup)
	window_size = int(args.wsize)
	quality_filter = 500
	snp_density = extract_vcf_data(vcf_file, quality_filter, '')
	mean_cov = extract_pileup_data(pileup_file)
	x, y, z = [], [], []
	fig = plt.figure(figsize=(15, 10))
	for element in snp_density:
		if element in mean_cov and mean_cov[element] < 200 and mean_cov[element] > 10 and snp_density[element] < 500:
			x.append(snp_density[element])
			y.append(mean_cov[element])
			z.append(element.split('_')[0])
	return x, y, z

def get_df2 (xy_dataset, scaf):
	DF2 = []
	count = 0
	while count <= (len(xy_dataset)-1):
		entry = xy_dataset[count]
		if entry[0].find(scaf) > -1:
			DF2.append([entry[0], entry[1], entry[2]])
		count = count + 1
	df2 = ''
	df2 = pd.DataFrame.from_dict(DF2).rename(columns={1:'SNPs2', 2:'coverage2'})		
	df2.replace(-np.inf, np.nan)
	df2.fillna(0)
	count = 0
	
	return df2

def var_v_cov_per_scaf(vcf, pileup, scaflist):
	x, y, z = var_v_cov(vcf, pileup)
	vcf_file = open(vcf)
	pileup_file = open(pileup)
	window_size = int(args.wsize)
	quality_filter = 500
	mean_cov = extract_pileup_data(pileup_file)
	snp_density = extract_vcf_data(vcf_file, quality_filter, '')
	x_scaf, y_scaf, z_scaf = [], [], []
	xy = np.vstack([x,y])
	xy_dataset = []
	for element in snp_density:	
		if element in mean_cov and mean_cov[element] < 200 and mean_cov[element] > 10 and snp_density[element] < 500:
			x_scaf.append(snp_density[element])
			y_scaf.append(mean_cov[element])
			xy_dataset.append([element[:element.find('_')], snp_density[element], mean_cov[element]])
	df = pd.DataFrame.from_dict(data=xy_dataset).rename(columns={0:'scaffolds', 1:'SNPs', 2:'coverage'})
	df.replace(-np.inf, np.nan)
	df.fillna(0)
	
	for scaf in scaflist:
		xy_copy = xy_dataset[:]
		df2 = ''
		df2 = df[df['scaffolds'].isin([scaf])]
		#print scaf+" mean and median SNPs/"+str(window_size), df2.SNPs2.mean(), df2.SNPs2.median()
		#print scaf+" mean and median coverage", df2.coverage2.mean(), df2.coverage2.median()
		fig = plt.figure(figsize=(30, 30))
		
		f, ax = plt.subplots(figsize=(8, 8))
		ax.set_aspect("equal")
		ax = sns.kdeplot(df.SNPs, df.coverage, cmap="Reds", shade=False)
		ax = sns.scatterplot(df2.SNPs, df2.coverage, cmap="Blues", marker='.')


		'''
		graph = sns.jointplot(x=df.SNPs, y=df.coverage, kind="kde", color='b')
		graph.x = df2.SNPs2
		graph.y = df2.coverage2
		graph.plot_joint(plt.scatter, marker='x', c='r', s=50)
		
		g = gaussian_kde(xy)(xy)
		plt.figure(figsize=(15,10))
		fig, ax = plt.subplots()
		x1 = pd.Series(x, name='SNPs in '+str(window_size)+" base pairs")
		x2 = pd.Series(y, name='Coverage')
		x_scaf1 = pd.Series(x_scaf, name='SNPs in '+str(window_size)+" base pairs")
		y_scaf2 = pd.Series(y_scaf, name='Coverage')
		
		ax = sns.jointplot(x_scaf1, y_scaf2, kind='scatter', color='g', linewidth=0.4, marker="+")
		ax = sns.jointplot(x1,x2, height=7, space=0, kind='kde')
		'''
		plt.savefig(args.output+'_var_v_cov.'+'ws'+str(window_size)+"_"+scaf+'.png')
		plt.close()
		
	


scaflist = set()
pileup = open(args.pileup)
	
for i in pileup:
	scaflist.add(i.split()[0])
pileup.seek(0)


var_v_cov_per_scaf(args.vcf, args.pileup, scaflist)
