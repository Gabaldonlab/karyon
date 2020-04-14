#!/bin/python

import sys
import scipy.stats
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import seaborn as sns

vcf_file = open(sys.argv[1])

value_list = []
binomial_list = []

for line in vcf_file:
	if line[0] == "#": continue
	else:
		chunk = line.split()[-1]
		if chunk.split(":")[0] == "0/1":
			values = chunk.split(':')[1].split(',')
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
plt.savefig('fair.png')
