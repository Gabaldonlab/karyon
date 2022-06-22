#!/bin/python
import sys, numpy, os.path, re
import argparse
from Bio import SeqIO


def soap_recipee(report, name, output, flags, job, loc):

	if output[-1] == "/":
		output = output[:-1]
	if os.path.isdir(output) == False:
		os.mkdir(output)
	locspp = output + "/" + name
	location=''

	###Dictionaries containing all the parameters for the libraries
	insert_size_dict = {}
	library_size_dict = {}
	type_dict = {}
	partner_dict = {}
	phred_dict = {}
	format_dict = {}
	compression_type_dict = {}
	deviation_dict = {}

	for line in open(report):
		chunk, a = line.split(), line.split()[0]
		insert_size_dict[a] = int(chunk[1])
		deviation_dict[a] = int(float(chunk[2])%1)
		library_size_dict[a] = int(chunk[3])
		phred_dict[a] = str(chunk[4])
		type_dict[a] = str(chunk[5])
		partner_dict[a] = chunk[6]
		format_dict[a] = chunk[7]
		compression_type_dict[a] = chunk[8]


	###Uses library info to write soapdenovo's config file###

	if os.path.exists(output+"/"+name+".soapdenovo_config.txt"):
		os.remove(output+"/"+name+".soapdenovo_config.txt")

	for fastqfile in insert_size_dict:
		soap_recipee=open(output+"/"+name+".soapdenovo_config.txt",'a')
		soap_recipee.write("[LIB]\n#average insert size\n")
		if type_dict[fastqfile] == '1':
			soap_recipee.write("avg_ins=" + str((insert_size_dict[fastqfile])*2)+"\n")
		elif type_dict[fastqfile] == '2': continue
		else:
			soap_recipee.write("avg_ins=" + str((insert_size_dict[fastqfile]))+"\n")
		
		soap_recipee.write("#if sequence needs to be reversed\n")
		if type_dict[fastqfile] == "pb":
			soap_recipee.write("reverse_seq=1\n")
		else:
			soap_recipee.write("reverse_seq=0\n")
		soap_recipee.write("#in which part(s) the reads are used\nasm_flags=3\n#use only first 100 bps of each read\n#rd_len_cutoff=100\n#in which order the reads are used while scaffolding\n")
		if type_dict[fastqfile] != '1' and type_dict[fastqfile] != '2':
			soap_recipee.write("rank=1\n")
		elif insert_size_dict[fastqfile] == '2': continue
		else:
			soap_recipee.write("rank=2\n")
		soap_recipee.write("#cutoff of pair number for a reliable connection (at least 3 for short insert size)\npair_num_cutoff=3\n#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\nmap_len=32\n#a pair of fastq file, read 1 file should always be followed by read 2 file\n")
		if type_dict[fastqfile] == '1':
			soap_recipee.write("q1="+fastqfile+"\n")
			soap_recipee.write("q2="+partner_dict[fastqfile]+"\n\n")
		elif insert_size_dict[fastqfile] == '2': continue
		else:	
			soap_recipee.write("q="+fastqfile+"\n\n")
	soap_recipee.close()


	if os.path.exists(output+"/"+name+".karyon.txt"):
		os.remove(output+"/"+name+".karyon.txt")

	#We produce a sh file that can be run with all the pipeline (soapdenovo)
	bash_job = job
	pairs = ''
	print (loc, "tomato")
	bash_job.write(loc+" all -s "+output+"/"+name+".soapdenovo_config.txt " + flags + " -o "+ locspp + " -p 12\n\n")


