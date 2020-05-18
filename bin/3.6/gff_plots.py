#!/bin/python
import numpy, math, sys, subprocess, os
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from decimal import Decimal
from scipy import stats
from scipy.stats import gaussian_kde
from Bio import SeqIO

vcf_file = pysam.VariantFile(sys.argv[1], 'r')
bam_file = pysam.AlignmentFile(sys.argv[2], 'rb')
fasta_file = sys.argv[3]
nQuire = "/root/src/karyon/src/dependencies/nQuire/nQuire"
kitchen = "/root/src/karyon/src/kitchen/"
window_size, step = 1000, 1000
counter = 20

def launch_nQuire(bam, nQuire, kitchen):
	BAMtemp = pysam.AlignmentFile(kitchen+"BAMtemp.bam", 'wb', template=bam_file)
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
	#print number_list
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

def window_walker(window_size, step, vcf_file, fasta_file, bam_file, nQuire, kitchen):
	window_stats = []
	prev_record_name = False
	fasta = SeqIO.parse(fasta_file, "fasta")
	for record in fasta:
		if record.name != prev_record_name and prev_record_name != False:
			nQuire_plot(window_stats)
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
				value = ref/alt
				if ref/alt == 0:
					value = float('nan')
				log_refalt_list.append(math.log(value,2))
			
			BAMtemp = bam_file.fetch(record.name, start, start + window_size)
			mean_cov = numpy.nanmean(bam_file.count_coverage(record.name, start, start + window_size, quality_threshold=0))

			stdev_cov = numpy.nanstd(bam_file.count_coverage(record.name, start, start + window_size))
			diplo_score, triplo_score, tetra_score = launch_nQuire(BAMtemp, nQuire, kitchen)
			R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB = ttest_ploidy(log_refalt_list)

			window_stats.append([window, start+window_size/2, diplo_score, triplo_score, tetra_score, R2_diploid, R2_triploidA, R2_triploidB, R2_tetraploidA, R2_tetraploidB, mean_cov, stdev_cov, mean_refaltcov_list])
			start = start + step
			
	return window_stats

def nQuire_plot(value_list):
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
		
	sns.set(style="darkgrid")

	fig, (p0,p1,p2,p3) = plt.subplots(nrows=4, sharex=True, figsize=(45,15))
	
	plt.subplot(4,1,1)
	plt.xlim(0,x[-1])
	plt.plot(x, y1, 'ro')
	plt.plot(x, y1, color='#aa0000', linestyle='--')
	plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
	plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

	plt.subplot(4,1,2)
	plt.xlim(0,x[-1])
	plt.plot(x, y2, 'go')
	plt.plot(x, y2, color='#00aa00', linestyle='--')
	plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
	plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

	plt.subplot(4,1,3)
	plt.xlim(0,x[-1])
	plt.plot(x, y3, 'bo')
	plt.plot(x, y3, color='#0000aa', linestyle='--')
	plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
	plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

	plt.subplot(4,1,4)
	plt.plot(x, mean_cov, color='grey', linestyle='-')


	#plt.subplot(5,1,5)
	#plt.plot(x, snp_den, 'm-')
	
	plt.subplot(5,1,4)
	plt.xlim(0,x[-1])
	plt.ylim(0,200)
	xy = np.vstack([pos_list,all_refalt_list])
	z = gaussian_kde(xy)(xy)
	plt.scatter(pos_list, all_refalt_list, c=z, s=30, edgecolor='')

	plt.savefig(name+".png")
	print(name+".png has been created")
	plt.clf()



value_list = window_walker(window_size, step, vcf_file, fasta_file, bam_file, nQuire, kitchen)
'''
x, y1, y2, y3, std_cov, mean_cov, snp_den = [], [], [], [], [], [], []
for i in value_list:
	if i[2] > 0:
		x.append(i[1]) 
		y1.append(i[2])
		y2.append(i[3])
		y3.append(i[4])
		std_cov.append(i[-2])
		mean_cov.append(i[-3])
		snp_den.append(i[-1])
	else: continue
		
sns.set(style="darkgrid")

fig, (p0,p1,p2,p3,p4) = plt.subplots(nrows=5, sharex=True, figsize=(45,15))
plt.subplot(5,1,1)
plt.plot(x, y1, 'r-')
plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

plt.subplot(5,1,2)
plt.plot(x, y2, 'g-')
plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

plt.subplot(5,1,3)
plt.plot(x, y3, 'b-')
plt.axhline(y=0.8, color='black', linestyle='-', linewidth=1)
plt.axhline(y=1, color='grey', linestyle='-', linewidth=1)

plt.subplot(5,1,4)
plt.plot(x, mean_cov, 'k-')
p3.errorbar(x, mean_cov, yerr=std_cov, fmt='b-')

plt.subplot(5,1,5)
plt.plot(x, snp_den, 'm-')

plt.savefig('za_plot.png')
'''

'''
def genome_ploidy_plot1 (nQuire_result):
	x = []
	y = []
	for i in nQuire_result:
		x.append(float(i.split()))
	fig, [ax1, ax2] = plt.subplots(2, 1, sharex=True)
	ax1.xcorr(x, y, usevlines=True, maxlags=50, normed=True, lw=2)
	ax1.grid(True)
	ax1.axhline(0, color='black', lw=2)

	ax2.acorr(x, usevlines=True, normed=True, maxlags=50, lw=2)
	ax2.grid(True)
	ax2.axhline(0, color='black', lw=2)

	plt.show()
'''

'''
def parse_pileup (pileup_file):
	coverage = []
	curr_pos = 0
	curr_scaffold = ''
	cov_dict = {}
	dev_dict = {}
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
	return cov_dict
'''	











