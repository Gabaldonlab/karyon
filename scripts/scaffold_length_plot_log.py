import sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

lenlist = []
fastainput = SeqIO.index(sys.argv[1], "fasta")

for i in fastainput:
	lenlist.append(len(fastainput.get_raw(i).decode()))

lenlist.sort()

N = len(lenlist)

ind = np.arange(N)
width = 0.35 

fig, ax = plt.subplots()
rects1 = ax.bar(ind, lenlist, width, color='royalblue')
ax.set_yscale('log')


plt.savefig('lenlog.png')
