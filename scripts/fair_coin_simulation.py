#!/bin/python

import sys
import scipy.stats
import numpy as np
import numpy.random
import matplotlib.pyplot as plt
import seaborn as sns

num = int(sys.argv[1])
cov = int(sys.argv[2])

data50 = []
data45raw = []
data33raw = []

for i in range(num):
	binomial50_value = (numpy.random.binomial(n=cov, p=0.5, size=None)/ float(cov))
	binomial52_value = (numpy.random.binomial(n=cov, p=0.55, size=None)/ float(cov))
	binomial48_value = (numpy.random.binomial(n=cov, p=0.45, size=None)/ float(cov))
	binomial33_value = (numpy.random.binomial(n=cov, p=0.33, size=None)/ float(cov))
	binomial66_value = (numpy.random.binomial(n=cov, p=0.66, size=None)/ float(cov))
	data50.append(binomial50_value)
	data45raw.append(binomial52_value)
	data45raw.append(binomial48_value)
	data33raw.append(binomial33_value)
	data33raw.append(binomial66_value)

data48 = numpy.random.choice(data45raw, size=(num/2), replace=False)
data33 = numpy.random.choice(data33raw, size=(num/2), replace=False)



bins = numpy.linspace(0, 1, 10000)

sns.distplot(data50, bins, hist=False, label='p=0.5', color="Maroon")
sns.distplot(data48, bins, hist=False, label='p=0.45', color="RoyalBlue")
sns.distplot(data33, bins, hist=False, label='p=0.33', color="Salmon")

plt.axvline(x=0.5, color='black', linestyle='-', linewidth=2)
plt.axvline(x=0.33, color='black', linestyle='--', linewidth=2)
plt.axvline(x=0.66, color='black', linestyle='--', linewidth=2)
plt.axvline(x=0.25, color='DimGray', linestyle='--', linewidth=1)
plt.axvline(x=0.75, color='DimGray', linestyle='--', linewidth=1)
plt.legend(loc='upper right')
plt.show()
