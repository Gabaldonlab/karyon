import sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt


lenlist = {}
fastainput = SeqIO.index(sys.argv[1], "fasta")
pileup_file = open(sys.argv[2])


def extract_pileup_data (pileup_file):
	x = []
	y = []
	coverage = []
	curr_pos = 0
	curr_scaffold = ''
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
	return x, y

x, y = extract_pileup_data(pileup_file)


plt.plot(x,y, '.')
ax = plt.subplot()
#ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')


plt.show()
plt.savefig('len_v_cov.png')





