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
from scipy.stats import gaussian_kde
import seaborn as sns
import pandas as pd
import scipy.stats
from decimal import Decimal
from scipy import stats

'''
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--fasta', required=True, help="fasta file used as input")
	parser.add_argument('-o', '--output', required=True, help="Output prefix")
	parser.add_argument('-v', '--vcf', required=True, help="VCF file used as input")
	parser.add_argument('-p', '--pileup', required=True, help="Mpileup file used as input")
	parser.add_argument('-b', '--bam', required=True, help="Bam file used as input")
	parser.add_argument('-l', '--library', required=True, help="Illumina library used for the KAT plot")
	parser.add_argument('-w', '--wsize', required=True, help="Window size for plotting")
	parser.add_argument('-c', '--counter', default=20, help="Number of scaffolds to analyze")
	parser.add_argument('-s', '--scafsize', default=False, help="Will ignore scaffolds with length below the given threshold")
	parser.add_argument('-S', '--scafmaxsize', default=False, help="Will ignore scaffolds with length above the given threshold")
	
args = parser.parse_args()
'''





def scaffold_len_lin(fasta, window_size, fastainput, output):
	lenlist = []
	for i in fastainput:
		lenlist.append(len(fastainput.get_raw(i).decode()))
	lenlist.sort()
	N = len(lenlist)
	ind = np.arange(N)
	width = 1 
	fig, ax = plt.subplots()
	ax.set_label('Scaffold length in nucleotides')
	rects1 = ax.bar(ind, lenlist, width, color='goldenrod')
	plt.axhline(y=window_size, color='dimgrey', linestyle='--', linewidth=2)
	plt.ylabel('Scaffold length')
	plt.savefig(output+'_lenlin.png')
	plt.close()

def scaffold_len_log(fasta, window_size, fastainput, output):
	lenlist = []
	for i in fastainput:
		lenlist.append(len(fastainput.get_raw(i).decode()))
	lenlist.sort()
	N = len(lenlist)
	ind = np.arange(N)
	width = 1
	fig, ax = plt.subplots()
	rects1 = ax.bar(ind, lenlist, width, color='royalblue')
	ax.set_yscale('log')
	plt.axhline(y=window_size, color='dimgrey', linestyle='--', linewidth=2)
	plt.ylabel('Scaffold length')
	plt.savefig(output+'_lenlog.png')
	plt.close()

def extract_vcf_data (vcf_file, quality_filter, window_size):
	snp_count = 0
	curr_pos = 0
	curr_scaffold = ''
	snp_dict = {}
	for line in vcf_file:
		if line[0] == "#": continue
		chunk = line.split("\t")
		#if args.scafsize != False:
			#if (lendict[chunk[0]]) <= scaf_size:continue
		#if args.scafmaxsize != False:
			#if (lendict[chunk[0]]) >= scaf_max_size: continue
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

def extract_pileup_data (pileup_file, window_size):
	coverage = []
	curr_pos = 0
	curr_scaffold = ''
	cov_dict = {}
	for line in pileup_file:
		if line[0] == "#": continue
		chunk = line.split()
		#if args.scafsize != False:
			#if (lendict[chunk[0]]) <= scaf_size:continue
		#if args.scafmaxsize != False:
			#if (lendict[chunk[0]]) >= scaf_max_size: continue
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
			coverage = [int(chunk[7][chunk[7].find('=')+1: chunk[7].find(';')])]
			#coverage = [ int(chunk[7].split(';')[0].split('=')[1]) ]
		else: 
			coverage.append(int(chunk[3]))
	pileup_file.seek(0)
	return cov_dict

def var_v_cov(vcf, pileup, window_size, output):
	vcf_file = open(vcf)
	pileup_file = open(pileup)
	quality_filter = 500
	snp_density = extract_vcf_data(vcf_file, quality_filter, window_size)
	mean_cov = extract_pileup_data(pileup_file, window_size)
	x, y = [], []
	fig = plt.figure(figsize=(15, 10))
	for element in snp_density:
		if element in mean_cov and mean_cov[element] < 200 and mean_cov[element] > 10 and snp_density[element] < 500:
			x.append(snp_density[element])
			y.append(mean_cov[element])
	xy = np.vstack([x,y])
	z = gaussian_kde(xy)(xy)
	plt.figure(figsize=(15,10))
	fig, ax = plt.subplots()
	x1 = pd.Series(x, name='SNPs in '+str(window_size)+" base pairs")
	x2 = pd.Series(y, name='Coverage')
	#ax.scatter(x, y, c=z, s=10, edgecolor='')
	ax = sns.jointplot(x1,x2,kind="kde", height=7, space=0)
	plt.savefig(output+'_var_v_cov.'+'ws'+str(window_size)+'.png')
	#plt.xlabel('SNPs in '+str(window_size)+"base pairs")
	#plt.ylabel('Coverage')
	plt.close()
	return mean_cov

def cov_plot(pileup, window_size, output):
	pileup_file = open(pileup)
	coverage = []
	curr_pos = 0
	curr_scaffold = ''
	cov_dict = {}
	fig = plt.figure(figsize=(15, 10))
	for line in pileup_file:
		if line[0] == "#": continue
		chunk = line.split()
		if chunk[0] != curr_scaffold:
			curr_scaffold = chunk[0]
			curr_pos = 0
			cov_dict[curr_scaffold+"_"+str(curr_pos)] = coverage
			coverage = []
		elif int(chunk[1]) >= curr_pos + window_size:
			cov_dict[curr_scaffold+"_"+str(curr_pos)] = coverage
			if (int(chunk[1])-curr_pos)/window_size > 1:
				for i in range(1, (int(chunk[1])-curr_pos)/window_size):
					cov_dict[curr_scaffold+"_"+str(curr_pos+(i*window_size))] = 0
			curr_pos = curr_pos + (window_size*((int(chunk[1])-curr_pos)/window_size))
			coverage = [int(chunk[3])]
		else: 
			coverage.append(int(chunk[3]))
	data = []
	meanlist = []
	for element in cov_dict:
		if numpy.mean(cov_dict[element]) < 200:
			data.append(cov_dict[element])
			meanlist.append(numpy.mean(cov_dict[element]))
		if len(data) >= 100:
			plt.boxplot(data[:100])
			break
	else:
		plt.boxplot(data)
	mean = numpy.mean(meanlist)
	plt.axhline(y=mean, color='dimgrey', linestyle='--', linewidth=2)
	plt.savefig(output+'_covplot.'+'ws'+str(window_size)+'.png')
	plt.close()
	pileup_file.close()

def fair_coin_global(vcf, window_size, output):
	vcf_file = open(vcf)
	value_list = []
	binomial_list = []
	fig = plt.figure(figsize=(20, 15))
	for line in vcf_file:
		if line[0] == "#": continue
		else:
			chunk = line.split()[-1]
			if chunk.split(":")[0] == "0/1":
				values = chunk.split(':')[1].split(',')
				if (int(values[0])+int(values[1])) == 0:
					value_list.append(float(values[1])/1)
					binomial_list.append(numpy.random.binomial(n=(int(values[0])+int(values[1])), p=0.5, size=None)/1)
				else:
					value_list.append(float(values[1])/(int(values[0])+int(values[1])))
					binomial_list.append(numpy.random.binomial(n=(int(values[0])+int(values[1])), p=0.5, size=None)/ (float(values[0])+int(values[1])))	
	#expected_freq = numpy.random.binomial(n=(sum(total_list)), p=0.5, size=None)
	#numpy.random.normal(loc=0.5, scale=np.std(value_list), size=len(value_list))

	print scipy.stats.chisquare(value_list, f_exp=binomial_list, ddof=0, axis=0)
	bins = numpy.linspace(0, 1, 100)
	sns.distplot(binomial_list, bins, label='exp')
	sns.distplot(value_list, bins, label='obs')
	plt.axvline(x=0.5, color='black', linestyle='-', linewidth=2)
	plt.axvline(x=0.33, color='black', linestyle='--', linewidth=2)
	plt.axvline(x=0.66, color='black', linestyle='--', linewidth=2)
	plt.axvline(x=0.25, color='DimGray', linestyle='--', linewidth=1)
	plt.axvline(x=0.75, color='DimGray', linestyle='--', linewidth=1)
	plt.legend(loc='upper right')
	plt.xlabel('Alternative/Reference SNP')
	plt.ylabel('Frequency')
	plt.savefig(output+'_faircoin_global.'+'ws'+str(window_size)+'.png')
	vcf_file.seek(0)
	plt.close()

def fair_coin_scaff(vcf, window_size, counter, output):
	vcf_file = open(vcf)
	scaffold_number = counter
	value_dict, alt_dict = {}, {}
	size_list = []
	binomial_list, binomial_altlist = [], []
	binomial_sublist, binomial_altsublist = [], []
	cutoff = 40
	fig = plt.figure(figsize=(15, 10))
	for line in vcf_file:
		if line[0] == "#": continue
		else:
			chunk = line.split()
		if chunk[0] not in value_dict:
			if len(value_dict) >= scaffold_number: break
			else:
				value_dict[chunk[0]] = []
				alt_dict[chunk[0]+"_alt"] = []
		if chunk[-1].split(":")[0] == "0/1":
			values = chunk[-1].split(':')[1].split(',')
			if int(values[0])+int(values[1]) >= cutoff:
				value_dict[chunk[0]].append(float(values[1])/(int(values[0])+int(values[1])))
				alt_dict[chunk[0]+"_alt"].append(float(values[0])/(int(values[0])+int(values[1])))
				binomial_value = numpy.random.binomial(n=(int(values[0])+int(values[1])), p=0.5, size=None)/ (float(values[0])+int(values[1]))
				binomial_list.append(binomial_value)

	for i in value_dict:
		size_list.append(len(value_dict[i]))
		sampling = int(np.mean(size_list))
		binomial_sublist = numpy.random.choice(binomial_list, size=sampling, replace=False)
		for i in binomial_sublist:
			binomial_altsublist.append(1-i)

		#expected_freq = numpy.random.binomial(n=(sum(total_list)), p=0.5, size=None)
		#numpy.random.normal(loc=0.5, scale=np.std(value_list), size=len(value_list))

		#print scipy.stats.chisquare(value_list, f_exp=binomial_list, ddof=0, axis=0)

		bins = numpy.linspace(0, 1, 10000)
	for value in value_dict:
		if len(value_dict[value]) > 0:
			sns.distplot(value_dict[value], bins, hist=False, color='RoyalBlue', norm_hist = True)
			#sns.distplot(alt_dict[value+"_alt"], bins, hist=False, color='DarkGoldenRod')
	
	sns.distplot(binomial_sublist, bins, hist=False, label='exp', color="Maroon")
	#sns.distplot(binomial_sublist, bins, hist=False, label='altexp', color="Maroon")
	plt.axvline(x=0.5, color='black', linestyle='-', linewidth=2)
	plt.axvline(x=0.33, color='black', linestyle='--', linewidth=2)
	plt.axvline(x=0.66, color='black', linestyle='--', linewidth=2)
	plt.axvline(x=0.25, color='DimGray', linestyle='--', linewidth=1)
	plt.axvline(x=0.75, color='DimGray', linestyle='--', linewidth=1)
	plt.legend(loc='upper right')
	plt.xlabel('Alternative/Reference SNP')
	plt.ylabel('Frequency')
	plt.savefig(output+'_faircoin_by_scaff.'+'ws'+str(window_size)+'.png')
	vcf_file.seek(0)
	plt.close()

def cov_v_len(pileup, fastainput, output):
	pileup_file = open(pileup)
	x = []
	y = []
	coverage = []
	curr_pos = 0
	curr_scaffold = ''
	fig = plt.figure(figsize=(15, 10))
	for line in pileup_file:
		if line[0] == "#": continue
		chunk = line.split()
		if chunk[0] != curr_scaffold:
			if curr_scaffold == '':
				curr_scaffold = chunk[0]
				coverage = []
				coverage.append(int(chunk[3]))
				continue
			else:
				curr_scaffold = chunk[0]
				y.append(np.mean(coverage))
				x.append(len(fastainput.get_raw(chunk[0]).decode()))
				coverage = []
				coverage.append(int(chunk[3]))
		else: 
			coverage.append(int(chunk[3]))
	plt.plot(x,y, '.')
	ax = plt.subplot()
	ax.set_xscale("log", nonposx='clip')
	ax.set_yscale("log", nonposy='clip')
	plt.xlabel('Scaffold length')
	plt.ylabel('Average coverage')
	plt.savefig(output+'_len_v_cov.png')
	pileup_file.seek(0)

def launch_nQuire(bam, nQuire, kitchen):
	BAMtemp = pysam.AlignmentFile(kitchen+"BAMtemp.bam", 'wb', template=bam)
	for i in bam:
		BAMtemp.write(i)
	BAMtemp.close()
	pysam.index(kitchen+"BAMtemp.bam")
	
	os.system(nQuire+" create -b "+kitchen+"BAMtemp.bam -o "+kitchen+"nQuire_temp")
	os.system(nQuire+" lrdmodel "+kitchen+"nQuire_temp.bin > "+kitchen+"nQuire_temp.report")
	nQuire_report = kitchen+"nQuire_temp.report"
		
	free_score, diplo_score, triplo_score, tetra_score = 0.0, 0.0, 0.0, 0.0
	for line in open(nQuire_report):
		if line.find("free") == 0: free_score = float(line.split()[1])
		elif line.find("dipl") == 0: diplo_score = float(line.split()[1])
		elif line.find("trip") == 0: triplo_score = float(line.split()[1])
		elif line.find("tetr") == 0: tetra_score = float(line.split()[1])
		else: continue
	os.remove(kitchen+"BAMtemp.bam")
	os.remove(kitchen+"BAMtemp.bam.bai")
	os.remove(kitchen+"nQuire_temp.bin")
	if free_score < 0.001:
		free_score, diplo_score, triplo_Score, tetra_Score = float('nan'), float('nan'), float('nan'), float('nan')
	return round(diplo_score/free_score, 3) , round(triplo_score/free_score, 3), round(tetra_score/free_score, 3)
	
def ttest_ploidy(number_list):
	y_dip, y_tripA, y_tripB, y_tetraA, y_tetraB = [], [], [], [], []
	for i in range(len(number_list)):
		y_dip.append(0)
		y_tripA.append(0.301)
		y_tripB.append(-0.301)
		y_tetraA.append(0.477)
		y_tetraB.append(-0.477)
	mean_number_list = numpy.nanmean(number_list)
	R2_diploid = (stats.ttest_1samp(number_list, math.log(1.0,2),nan_policy='omit'))
	R2_triploidA, R2_triploidB  = (stats.ttest_1samp(number_list, math.log(2.0,2),nan_policy='omit')), (stats.ttest_1samp(number_list, math.log(2.0,2)*-1,nan_policy='omit'))
	R2_tetraploidA, R2_tetraploidB  = (stats.ttest_1samp(number_list, math.log(3.0,2),nan_policy='omit')), (stats.ttest_1samp(number_list, math.log(1/3.0,2),nan_policy='omit'))
	return R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB#, mean_number_list, numpy.nanstd(number_list)

def window_walker(window_size, step, vcf_file, fasta_file, bam_file, nQuire, kitchen, newpath, counter):
	window_stats = []
	prev_record_name = False
	fasta = SeqIO.parse(fasta_file, "fasta")
	for record in fasta:
		if len(os.listdir(newpath+"nQuireplots_ws"+str(window_size))) >= counter: break
		if record.name != prev_record_name and prev_record_name != False:
			nQuire_plot(window_stats, window_size, newpath)
			window_stats = []
		if prev_record_name == False:
			prev_record_name = record.name
		start, end = 0, len(record)
		log_refalt_list = []
		while start + step <= end:
			window = record.name+":"+str(start)+":"+str(start+window_size)
			VCF = vcf_file.fetch(record.name, start, start + window_size) 
			mean_refaltcov_list = []
			for i in VCF:
				refalt = (str(i).split()[-1].split(':')[1].split(","))
				snp_pos = int(str(i).split()[1])
				ref, alt = float(refalt[0]), float(refalt[1])
				totalcov = numpy.nanmean(bam_file.count_coverage(record.name, snp_pos, snp_pos+1, quality_threshold=0))
				mean_refaltcov_list.append(totalcov)
				if alt == 0 or ref == 0:
					value = float('nan')
				else:
					value = alt/(float(ref)+0.01)
				log_refalt_list.append(math.log(value,2))
			
			BAMtemp = bam_file.fetch(record.name, start, start + window_size)
			mean_cov = numpy.nanmean(bam_file.count_coverage(record.name, start, start + window_size, quality_threshold=0))

			stdev_cov = numpy.nanstd(bam_file.count_coverage(record.name, start, start + window_size))
			diplo_score, triplo_score, tetra_score = launch_nQuire(BAMtemp, nQuire, kitchen)
			R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB = ttest_ploidy(log_refalt_list)

			window_stats.append([window, start+window_size/2, diplo_score, triplo_score, tetra_score, R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB, mean_cov, stdev_cov, mean_refaltcov_list])
			start = start + step
			vcf_file.seek(0)
	return window_stats

def nQuire_plot(value_list, window_size, newpath):
	name, x, y1, y2, y3, std_cov, mean_cov, snp_den = '', [], [], [], [], [], [], []
	all_refalt_list, pos_list = [], []
	for i in value_list:
		if i[2] > 0:
			name = i[0].split(":")[0]
			x.append(i[1]) 
			if i[2] > 1: y1.append(1.0)
			else: y1.append(i[2])
			if i[3] > 1: y2.append(1.0)
			else: y2.append(i[3])
			if i[4] > 1: y3.append(1.0)
			else: y3.append(i[4])
			std_cov.append(i[-1])
			mean_cov.append(i[-3])
			snp_den.append(i[-2])
			for e in i[-1]:
				pos_list.append(i[1])
				all_refalt_list.append(e)
		else: continue
	if len(x) > 0:
		sns.set(style="darkgrid")
	
		fig, (p0,p1,p2,p3,p4) = plt.subplots(nrows=5, sharex=True, figsize=(45,15))
		plt.subplot(5,1,1)
		plt.xlim(0,x[-1])
		plt.plot(x, y1, 'ro')
		plt.plot(x, y1, color='#aa0000', linestyle='--')
		plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
		plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

		plt.subplot(5,1,2)
		plt.xlim(0,x[-1])
		plt.plot(x, y2, 'go')
		plt.plot(x, y2, color='#00aa00', linestyle='--')
		plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
		plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

		plt.subplot(5,1,3)
		plt.xlim(0,x[-1])
		plt.plot(x, y3, 'bo')
		plt.plot(x, y3, color='#0000aa', linestyle='--')
		plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
		plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

		plt.subplot(5,1,4)
		plt.plot(x, mean_cov, color='grey', linestyle='-')
		
		plt.subplot(5,1,5)
		plt.xlim(0,x[-1])
		plt.ylim(0,200)
		xy = np.vstack([pos_list,all_refalt_list])
		z = gaussian_kde(xy)(xy)
		plt.scatter(pos_list, all_refalt_list, c=z, s=30, edgecolor='')
		plt.savefig(newpath+"nQuireplots_ws"+str(window_size)+"/"+name+".png")
		print newpath+"nQuireplots_ws"+str(window_size)+"/"+name+".png has been created"
		plt.clf()
	

def katplot(fasta, library, KAT, out):
	# os.system(KAT+" comp -o "+out+" "+library+" "+fasta+" > "+out+".katreport")
	cmd = KAT+" comp -o "+out+" "+library+" "+fasta+" > "+out+".katreport"
	returned_value = subprocess.call(cmd, shell=True)  # returns the exit code in unix
	print '###############'
	print ('KAT:', returned_value)
	print '###############'


def allplots(window_size, vcf, fasta_file, bam, mpileup, library, nQuire, KAT, kitchen, newpath, counter, kitchenID, out_name):
	if out_name==False:
		outname = ''
	newpath = newpath+"/"+out_name
	# os.system("bgzip -c "+ vcf+ " > " + vcf + ".gz")
	# os.system("tabix -p vcf "+ vcf+".gz")
	vcf_file = open(vcf+".gz", 'r')
	bam_file = pysam.AlignmentFile(bam, 'rb')
	kitchen = kitchen+kitchenID

	lendict = {}
	fastainput = SeqIO.index(fasta_file, "fasta")
	for i in fastainput:
		lendict[i] = len(fastainput.get_raw(i).decode())
	step = window_size/2
	# VCF = pysam.VariantFile(vcf+".gz", 'r')
	if newpath.find("/") > -1:
		originalpath=os.getcwd()
	os.makedirs(newpath+"nQuireplots_ws"+str(window_size))
	
	scaffold_len_lin(fasta_file, window_size, fastainput, newpath)
	scaffold_len_log(fasta_file, window_size, fastainput, newpath)
	var_v_cov(vcf, mpileup, window_size, newpath)
	cov_plot(mpileup, window_size, newpath)
	fair_coin_global(vcf, window_size, newpath)
	fair_coin_scaff(vcf, window_size, counter, newpath)
	cov_v_len(mpileup, fastainput, newpath)
	katplot(fasta_file, library, KAT, newpath)
	# window_walker(window_size, step, VCF, fasta_file, bam, nQuire, kitchen, newpath, counter)

#newpath = args.output[:args.output.rfind("/")]+"/"
#allplots(window_size, vcf_file, fasta_file, bam_file, mpileup, library, nQuire, KAT, kitchen, newpath, counter, kitchenID)



	



	
	
	
