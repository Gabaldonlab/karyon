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
width = 1 

fig, ax = plt.subplots()
ax.set_label('Scaffold length in nucleotides')
rects1 = ax.bar(ind, lenlist, width, color='goldenrod')


plt.savefig('lenlin.png')
