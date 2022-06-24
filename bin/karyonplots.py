#!/bin/python
import sys, os, re, subprocess, math
import psutil
from pysam import pysam
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

def extract_vcf_data (vcf_file, quality_filter, window_size, lendict, scafminsize, scafmaxsize):
	snp_count = 0
	curr_pos = 0
	curr_scaffold = ''
	snp_dict = {}
	for line in vcf_file:
		if line[0] == "#": continue
		chunk = line.split("\t")
		if scafminsize != False and lendict[chunk[0]] <= scafminsize: continue
		elif scafmaxsize != False and lendict[chunk[0]] >= scafminsize: continue
		else:
			if float(chunk[5]) < quality_filter: continue
			if chunk[0] != curr_scaffold:
				curr_scaffold = chunk[0]
				curr_pos = 0
				snp_dict[curr_scaffold+"_"+str(curr_pos)] = snp_count
				snp_count = 1
			elif int(chunk[1]) >= curr_pos + window_size:
				snp_dict[curr_scaffold+"_"+str(curr_pos)] = snp_count
				if (int(chunk[1])-curr_pos)/window_size > 1:
					for i in range(1, int((int(chunk[1])-curr_pos)//window_size)):
						snp_dict[curr_scaffold+"_"+str(curr_pos+(i*window_size))] = 0
				curr_pos = curr_pos + (window_size*((int(chunk[1])-curr_pos)/window_size))
				snp_count = 1
			else:
				snp_count = snp_count + 1
	vcf_file.seek(0)
	return snp_dict

def extract_pileup_data (pileup_file, window_size, lendict, scafminsize, scafmaxsize):
	coverage = []
	curr_pos = 0
	curr_scaffold = ''
	cov_dict = {}
	for line in pileup_file:
		if scafminsize != False and lendict[chunk[0]] <= scafminsize: continue
		elif scafmaxsize != False and lendict[chunk[0]] >= scafminsize: continue
		else:
			if line[0] == "#": continue
			chunk = line.split()
			if chunk[0] != curr_scaffold:
				curr_scaffold = chunk[0]
				curr_pos = 0
				cov_dict[curr_scaffold+"_"+str(curr_pos)] = numpy.mean(coverage)
				coverage = []
			elif int(chunk[1]) >= curr_pos + window_size:
				cov_dict[curr_scaffold+"_"+str(curr_pos)] = numpy.mean(coverage)
				if (int(chunk[1])-curr_pos)/window_size > 1:
					for i in range(1, int((int(chunk[1])-curr_pos)//window_size)):
						cov_dict[curr_scaffold+"_"+str(curr_pos+(i*window_size))] = 0
				curr_pos = curr_pos + (window_size*((int(chunk[1])-curr_pos)/window_size))
				res = max(line, key = len)
				coverage.append(int(chunk[3]))
			else:
				res = max(line, key = len)
				coverage.append(int(chunk[3]))
	pileup_file.seek(0)
	return cov_dict

def var_v_cov(vcf, pileup, window_size, output, lendict, scafminsize, scafmaxsize):
	vcf_file = open(vcf)
	pileup_file = open(pileup)
	quality_filter = 500
	snp_density = extract_vcf_data(vcf_file, quality_filter, window_size, lendict, scafminsize, scafmaxsize)
	mean_cov = extract_pileup_data(pileup_file, window_size, lendict, scafminsize, scafmaxsize)
	x, y = [], []
	fig = plt.figure(figsize=(15, 10))
	for element in snp_density:
		if element in mean_cov and mean_cov[element] < 200 and mean_cov[element] > 10 and snp_density[element] < 100:
			x.append(snp_density[element])
			y.append(mean_cov[element])
	xy = np.vstack([x,y])
	plt.figure(figsize=(15,10))
	fig, ax = plt.subplots()
	x1 = pd.Series(x, name='SNPs in '+str(window_size)+" base pairs")
	x2 = pd.Series(y, name='Coverage')
	ax = sns.jointplot(x1,x2,kind="kde", height=7, space=0, legend=True, marginal_ticks=True, fill=True, bw_adjust=.5)
	plt.savefig(output+'_var_v_cov.'+'ws'+str(window_size)+'.png')
	plt.xlabel('SNPs in '+str(window_size)+"base pairs")
	plt.ylabel('Coverage')
	plt.close()
	return mean_cov

def cov_plot(pileup, window_size, output, lendict, scafminsize, scafmaxsize):
	pileup_file = open(pileup)
	coverage = []
	curr_pos = 0
	curr_scaffold = ''
	cov_dict = {}
	fig = plt.figure(figsize=(15, 10))
	for line in pileup_file:
		if line[0] == "#": continue
		chunk = line.split()
		if scafminsize != False and lendict[chunk[0]] <= scafminsize: continue
		elif scafmaxsize != False and lendict[chunk[0]] >= scafminsize: continue
		else:
			if chunk[0] != curr_scaffold:
				curr_scaffold = chunk[0]
				curr_pos = 0
				cov_dict[curr_scaffold+"_"+str(curr_pos)] = coverage
				coverage = []
			elif int(chunk[1]) >= curr_pos + window_size:
				cov_dict[curr_scaffold+"_"+str(curr_pos)] = coverage
				if (int(chunk[1])-curr_pos)/window_size > 1:
					for i in range(1, int((int(chunk[1])-curr_pos)//window_size)):
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

def fair_coin_global(vcf, window_size, output, lendict, scafminsize, scafmaxsize):
	vcf_file = open(vcf)
	value_list = []
	binomial_list = []
	fig = plt.figure(figsize=(20, 15))
	for line in vcf_file:
		if line[0] == "#": continue
		else:
			chunk = line.split()[-1]
			if scafminsize != False and lendict[chunk[0]] <= scafminsize: continue
			elif scafmaxsize != False and lendict[chunk[0]] >= scafminsize: continue
			else:
				if chunk.split(":")[0] == "0/1":
					values = chunk.split(':')[1].split(',')
					if (int(values[0])+int(values[1])) == 0:
						value_list.append(float(values[1])/1)
						binomial_list.append(numpy.random.binomial(n=(int(values[0])+int(values[1])), p=0.5, size=None)/1)
					else:
						value_list.append(float(values[1])/(int(values[0])+int(values[1])))
						binomial_list.append(numpy.random.binomial(n=(int(values[0])+int(values[1])), p=0.5, size=None)/ (float(values[0])+int(values[1])))	
	#print (scipy.stats.chisquare(value_list, f_exp=binomial_list, ddof=0, axis=0))
	bins = numpy.linspace(0, 1, 100)
	sns.displot(binomial_list, bins, label='exp')
	sns.displot(value_list, bins, label='obs')
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

def fair_coin_scaff(vcf, window_size, counter, output, lendict, scafminsize, scafmaxsize):
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
		if scafminsize != False and lendict[chunk[0]] <= scafminsize: continue
		elif scafmaxsize != False and lendict[chunk[0]] >= scafminsize: continue
		else:	
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
		bins = numpy.linspace(0, 1, 10000)
	for value in value_dict:
		if len(value_dict[value]) > 0:
			sns.displot(value_dict[value], bins, hist=False, color='RoyalBlue', norm_hist = True)	
	sns.displot(binomial_sublist, bins, hist=False, label='exp', color="Maroon")
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
	ax.set_xscale("log", nonpositive='clip')
	ax.set_yscale("log", nonpositive='clip')
	plt.xlabel('Scaffold length')
	plt.ylabel('Average coverage')
	plt.savefig(output+'_len_v_cov.png')
	pileup_file.seek(0)

def launch_nQuire(BAMtemp, nQuire, kitchen):
	
	os.system(nQuire+" create -b "+kitchen+"BAMtemp.bam -o "+kitchen+"nQuire_temp")
	os.system(nQuire+" lrdmodel "+kitchen+"nQuire_temp.bin > "+kitchen+"nQuire_temp.report")
	nQuire_report = kitchen+"nQuire_temp.report"
		
	free_score, diplo_score, triplo_score, tetra_score = 0.0, 0.0, 0.0, 0.0
	for line in open(nQuire_report):
		if line.find("file") == 0: continue
		else:
			free_score = float(line.split()[1])
			diplo_score = float(line.split()[2])
			triplo_score = float(line.split()[3])
			tetra_score = float(line.split()[4])

	if free_score < 0.001:
		free_score, diplo_score, triplo_Score, tetra_Score = float('nan'), float('nan'), float('nan'), float('nan')
	return round(diplo_score/free_score, 3) , round(triplo_score/free_score, 3), round(tetra_score/free_score, 3), free_score
	
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
	return R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB

def window_walker(window_size, step, vcf, fasta_file, bam, nQuire, kitchen, newpath, counter, lendict, scafminsize, scafmaxsize, no_plot):
	N, dfbatch = 0, []
	step = int(step)
	bam_df = pd.DataFrame([x.split('\t') for x in pysam.depth(bam).split('\n')])
	vcf_file = pysam.VariantFile(vcf+".gz", 'r')
	bam_file = pysam.AlignmentFile(bam, 'rb')
	vcfset = set()
	if newpath.find("/") > -1:
		originalpath=os.getcwd()
	if os.path.isdir(newpath+"nQuireplots_ws"+str(window_size)):
		os.rmdir(newpath+"nQuireplots_ws"+str(window_size))
	os.makedirs(newpath+"nQuireplots_ws"+str(window_size))
	window_stats = []
	prev_record_name = False
	fasta = SeqIO.parse(fasta_file, "fasta")
	for record in fasta:
		if len(os.listdir(newpath+"nQuireplots_ws"+str(window_size))) >= counter: break
		if record.name != prev_record_name and prev_record_name != False:
			n = 0
			cov_list = [[],[]]
			while n+step <= end:
				cov_list[0].append(n+step/2)
				mask = bam_df[bam_df[0] == record.name][2][n:n+step].astype(int)
				cov_list[1].append(mask.mean())
				#a = pysam.depth("-aa", "-r", record.name+":"+str(n)+"-"+str(n+step), bam).split()
				#if len(a) == 0:
				#	cov_list[1].append(0.0)
				#else:
				#	cov_list[1].append(float(pysam.depth("-aa", "-r", record.name+":"+str(n+1)+"-"+str(n+step), bam).split()[-1]))
				n = n + step

			if no_plot == False:
				nQuire_plot(window_stats, window_size, newpath, cov_list[0], cov_list[1], lendict, scafminsize, scafmaxsize)
			window_stats = []
		if prev_record_name == False:
			prev_record_name = record.name
		start, end = 0, len(record)
		log_refalt_list = []
		while start + step <= end:
			window = record.name+":"+str(start)+":"+str(start+window_size)
			VCF = vcf_file.fetch(contig=record.name, start=start, stop=start + window_size) 
			mean_refaltcov_list = []
			for i in VCF:
				refalt = (str(i).split()[-1].split(':')[1].split(","))
				snp_pos = int(str(i).split()[1])
				ref, alt = float(refalt[0]), float(refalt[1])
				A = bam_file.count_coverage(record.name, snp_pos, snp_pos+1, quality_threshold=0)
				totalcov = sum([A[0][0], A[1][0], A[2][0], A[3][0]])
				mean_refaltcov_list.append(totalcov)
				if alt == 0 or ref == 0:
					value = float('nan')
				else:
					value = alt/(float(ref)+0.01)
				log_refalt_list.append(math.log(value,2))
			BAMfetch = bam_file.fetch(record.name, start, start + window_size)
			BAMtemp = pysam.AlignmentFile(kitchen+"BAMtemp.bam", 'wb', template=bam_file)
			for i in BAMfetch:
				BAMtemp.write(i)
			BAMtemp.close()
			pysam.index(kitchen+"BAMtemp.bam")
			mean_cov = numpy.nanmean(bam_file.count_coverage(record.name, start, start + window_size, quality_threshold=0))

			stdev_cov = numpy.nanstd(bam_file.count_coverage(record.name, start, start + window_size))
			diplo_score, triplo_score, tetra_score, free_score = launch_nQuire(BAMtemp, nQuire, kitchen)
			R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB = ttest_ploidy(log_refalt_list)
			df = pd.DataFrame(columns=["contig", "start", "end", "diplo_score", 
			"triplo_score", "tetra_score", "free_score", "normalized_diplo_score", "normalized_triplo_score", "normalized_tetra_score", "mean_cov", 
			"st_dev_cov", "SNPs", "ploidy"])
			window_stats.append([window, start+(window_size/2), diplo_score, triplo_score, tetra_score, R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB, mean_cov, stdev_cov, mean_refaltcov_list])
			dfrow = pd.DataFrame({"contig":record.name, "start": start, "end": start+window_size, "diplo_score": diplo_score*free_score, "triplo_score": triplo_score*free_score, "tetra_score": tetra_score*free_score,
			"free_score": free_score, "normalized_diplo_score":diplo_score, "normalized_triplo_score": triplo_score, "normalized_tetra_score": tetra_score,"mean_cov": mean_cov, "st_dev_cov": stdev_cov, "SNPs": len(mean_refaltcov_list), "ploidy":0}, index=[N])
			start = start + step
			dfbatch.append(dfrow)
			N = N+1
		vcf_file.seek(0)
	df = pd.concat(dfbatch, ignore_index=True)
	return(df)

def nQuire_plot(value_list, window_size, newpath, xcov, ycov, lendict, scafminsize, scafmaxsize):
	name, x, y1, y2, y3, std_cov, mean_cov, snp_den = '', [], [], [], [], [], [], []
	all_refalt_list, pos_list = [], []
	for i in value_list:
		if i[2] > 0:
			name = i[0].split(":")[0].replace("|", "_")
			if scafminsize != False and lendict[name] <= scafminsize: continue
			elif scafmaxsize != False and lendict[name] >= scafminsize: continue
			else:
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

	if len(x) > 0:
		if scafminsize != False and lendict[name] <= scafminsize:
			pass
		elif scafmaxsize != False and lendict[name] >= scafminsize:
			pass
		else:
			sns.set(style="darkgrid")
		
			#Diploid score
			fig, (p0,p1,p2,p3,p4) = plt.subplots(nrows=5, sharex=True, figsize=(45,15))
			plt.title(name)
			plt.subplot(5,1,1)
			plt.xlim(0,x[-1])
			plt.plot(x, y1, 'ro')
			plt.plot(x, y1, color=(0.7,0,0), linestyle='--')
			plt.axhline(y=0.8, color=(0,0,0), linestyle='-', linewidth=1)
			plt.axhline(y=1, color=(0.3,0.3,0.3), linestyle='-', linewidth=1)
			plt.ylabel('Diploid score')
		
			#Triploid score
			plt.subplot(5,1,2)
			plt.xlim(0,x[-1])
			plt.plot(x, y2, 'go')
			plt.plot(x, y2, color=(0, 0.7, 0), linestyle='--')
			plt.axhline(y=0.8, color=(0,0,0), linestyle='-', linewidth=1)
			plt.axhline(y=1, color=(0.3, 0.3, 0.3), linestyle='-', linewidth=1)
			plt.ylabel('Triploid score')
		
			#Tetraploid score
			plt.subplot(5,1,3)
			plt.xlim(0,x[-1])
			plt.plot(x, y3, 'bo')
			plt.plot(x, y3, color=(0, 0, 0.7), linestyle='--')
			plt.axhline(y=0.8, color=(0,0,0), linestyle='-', linewidth=1)
			plt.axhline(y=1, color=(0.3, 0.3, 0.3), linestyle='-', linewidth=1)
			plt.ylabel('Tetraploid score')

			#Coverage
			plt.subplot(5,1,4)
			plt.xlim(0,x[-1])
			plt.plot(xcov, ycov, color=(0.3, 0.3, 0.3), linestyle='-')
			plt.ylabel('Depth of coverage')
		
			#SNP location
			plt.subplot(5,1,5)
			plt.xlim(0,x[-1])
			plt.ylim(0,200)
			xy = np.vstack([pos_list,all_refalt_list])
			plt.ylabel('SNP location and coverage')
			if len(xy[0]) > 3:
				if len(set(xy[0])) == 1:
					xy[0][0] = xy[0][0]+10
					xy[0][-1] = xy[0][-1]-10
					z = gaussian_kde(xy, bw_method=0.001)(xy)
				else:
					z = gaussian_kde(xy)(xy)
				plt.scatter(pos_list, all_refalt_list, c=z, s=30, edgecolor=(0.4, 0.4, 0.4))
			plt.savefig(newpath+"nQuireplots_ws"+str(window_size)+"/"+name+".png")
			print (newpath+"nQuireplots_ws"+str(window_size)+"/"+name+".png has been created")
			plt.clf()	

def katplot(fasta, library, KAT, out):
	cmd = "conda run -n kat_env kat comp -o "+out+" "+library+" "+fasta+" > "+out[:-1]+".katreport"
	returned_value = subprocess.call(cmd, shell=True)  # returns the exit code in unix
	print ('###############')
	print ('KAT:', returned_value)
	print (cmd)
	print ('###############')

def allplots(window_size, vcf, fasta_file, bam, mpileup, library, nQuire, KAT, kitchen, newpath, counter, kitchenID, out_name, scafminsize, scafmaxsize, no_plot):
	if out_name==False:
		outname = ''
	if newpath[-1] != "/":
		newpath = newpath + "/"
	newpath = newpath+out_name
	os.system("bgzip -c "+ vcf+ " > " + vcf + ".gz")
	os.system("tabix -p vcf "+ vcf+".gz")
	vcf_file = open(vcf, 'r')
	bam_file = pysam.AlignmentFile(bam, 'rb')
	if os.path.exists(kitchen) == False:
		if os.path.exists(kitchen[:kitchen.rfind("/")]) == False:
			os.makedirs(kitchen[:kitchen.rfind("/")])
			if os.path.exists(kitchen) == False:
				os.makedirs(kitchen)
	lendict = {}
	fastainput = SeqIO.index(fasta_file, "fasta")
	for i in fastainput:
		lendict[i] = len(fastainput.get_raw(i).decode())
	step = int(window_size/2)
	df = window_walker(window_size, step, vcf, fasta_file, bam, nQuire, kitchen, newpath, counter, lendict, scafminsize, scafmaxsize, no_plot)
	if no_plot == False:
		scaffold_len_lin(fasta_file, window_size, fastainput, newpath)
		scaffold_len_log(fasta_file, window_size, fastainput, newpath)
		var_v_cov(vcf, mpileup, window_size, newpath, lendict, scafminsize, scafmaxsize)
		cov_plot(mpileup, window_size, newpath, lendict, scafminsize, scafmaxsize)
		fair_coin_global(vcf, window_size, newpath, lendict, scafminsize, scafmaxsize)
		fair_coin_scaff(vcf, window_size, counter, newpath, lendict, scafminsize, scafmaxsize)
		cov_v_len(mpileup, fastainput, newpath)
		katplot(fasta_file, library, KAT, newpath)
	return(df)
	
	



	
	
	
