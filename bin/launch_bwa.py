# BWA for map in chloroplast
# Irene Julca 29/06/15

import argparse, os, os.path
from os import listdir
from os.path import isfile, join


parser = argparse.ArgumentParser(description="BWA mapping for one fastq file (sam and bam files).")
parser.add_argument("-r", "--ref", dest="ref", required=True, help="genome of reference")
parser.add_argument("-f1", "--f1path", dest="f1path", required=True, help="fastq_1 file")
parser.add_argument("-f2", "--f2path", dest="f2path", default=False, help="fastq_2 file")
parser.add_argument("-n", "--name", dest="outname", required=True, help="name of the output")
args = parser.parse_args()

fasq1 = str(args.f1path)
fasq2 = str(args.f2path)
reference = str(args.ref)
name = str(args.outname)

if args.f2path == False:
	comando1 = 'bwa mem -R "@RG\\tID:test\\tSM:bar" -t 8 ' + reference + ' ' + fasq1 + " " + ">" + name + ".sam"
else:
	comando1 = 'bwa mem -R "@RG\\tID:test\\tSM:bar" -t 8 ' + reference + ' ' + fasq1 + " " + fasq2 + " >" + name + ".sam"
os.system(comando1)

comando2 = "samtools view -Sb " + name + ".sam >" + name + ".bam"
os.system(comando2)

comando3 = "samtools sort " + name + ".bam -T " + name + " >" + name + ".sorted.bam"
os.system(comando3)

comando4 = "samtools index " + name + ".sorted.bam"
os.system(comando4)

print("...")
print("...")	 
		
