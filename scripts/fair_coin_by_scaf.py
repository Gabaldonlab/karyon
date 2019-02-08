#!/bin/python

import sys
import scipy.stats
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import seaborn as sns

vcf_file = open(sys.argv[1])
scaffold_number = int(sys.argv[2])

value_dict, alt_dict = {}, {}
size_list = []
binomial_list, binomial_altlist = [], []
binomial_sublist, binomial_altsublist = [], []
cutoff = 40

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
print sampling
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
plt.savefig('fair2.png')
