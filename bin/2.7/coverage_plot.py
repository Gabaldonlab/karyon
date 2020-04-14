#!/bin/python

import sys, numpy
from Bio import SeqIO
import matplotlib.pyplot as plt
from operator import itemgetter

pileup_file = open(sys.argv[1])
window_size = int(sys.argv[2])
fastafile = SeqIO.index(sys.argv[3], "fasta")

'''lendict = {}
for i in fastafile:
	chunk = i.split("size")
	newnum = int(chunk[1])/window_size
	newnum = newnum * window_size
	newname = i + "_" + str(newnum)
	lendict[newname] = len(fastafile.get_raw(i).decode())

print lendict
'''
quality_filter = 0

def extract_pileup_data (pileup_file):
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
	
fig = plt.figure(figsize=(15, 10))
mean_cov = extract_pileup_data(pileup_file)

data = []
meanlist = []

for element in mean_cov:
	if numpy.mean(mean_cov[element]) < 1000:
		data.append([mean_cov[element], numpy.mean(mean_cov[element])])
		meanlist.append([numpy.mean(mean_cov[element]), numpy.mean(mean_cov[element])])
		data = sorted(data, key=itemgetter(1), reverse=False)
		meanlist = sorted(meanlist, key=itemgetter(1), reverse=False)
		sortdata, sortmeanlist = [], []
		for i in data:
			sortdata.append(i[0])
		for i in meanlist:
			sortdata.append(i[0])
if len(data) >= 100:
	plt.boxplot(sortdata[:100])
else:
	plt.boxplot(sortdata)
mean = numpy.mean(sortmeanlist)
plt.axhline(y=mean, color='dimgrey', linestyle='--', linewidth=2)
plt.savefig('covplot.png')





