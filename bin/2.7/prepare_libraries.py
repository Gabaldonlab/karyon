#!/bin/python
import sys, numpy, os.path, re, os
import argparse
from Bio import SeqIO
import gzip, bz2, tarfile

def remove_false_files(filelist): #We remove empty files that otherwise would make everything crash
	clean_list = []
	for i in filelist: 
		if os.path.getsize(i) < 10000: #I probably need a better way to judge it
			continue
		else:
			clean_list.append(i)
	return clean_list

def set_sampling(fastqfile, sample_size):
	sampling = []
	fastafile = SeqIO.parse(fastqfile, "fasta")
	for i in fastafile:
		sampling.append(len(i))
		if len(sampling) >= sample_size:
			break
	return sampling

def get_mean_read_len(fastqfile, sample_size, compressed_dict):

	sampling = []
	if (compressed_dict[fastqfile] != "no-compression" and fastqfile[:fastqfile.rfind(".")+1][-3:] == "fa") or (compressed_dict[fastqfile] != "no-compression" and  fastqfile[:fastqfile.rfind(".")+1][-6:] == ".fasta"):
		sampling = set_sampling(fastqfile, sample_size)
	elif fastqfile[-3:] == "fa" or fastqfile[-6:] == ".fasta" :
		sampling = set_sampling(fastqfile, sample_size)
	else:
		switch = False
		if compressed_dict[fastqfile] == "gzip":
			open_fastqfile = gzip.open(fastqfile, 'r')
		elif compressed_dict[fastqfile] == "bz2":
			open_fastqfile = bz2.BZ2file(fastqfile, "r")
		elif compressed_dict[fastqfile] == "tar":
			open_fastqfile = tarfile.open(fastqfile, "r")
		else:
			open_fastqfile = open(fastqfile)
		for line in open_fastqfile:
			if switch == True:
				sampling.append(len(line)-1)
				switch = False
				if len(sampling) >= sample_size:
					break
			elif line[0] == "+":
				switch = True
			else: continue
	return [int(numpy.mean(sampling)), numpy.std(sampling)]

def compression_parse(fastq):
	compressed_dict = {}
	for i in fastq:
		if i[i.rfind("gz"):] == "gz" or i[i.rfind("gzip"):] == "gzip":
			compressed_dict[i] = "gzip"
		elif i[i.rfind("bzip2"):] == "bz2" or i[i.rfind("bzip2"):] == "bzip2" > -1:
			compressed_dict[i] = "bz2"
		elif i[i.rfind("tar"):] == "tar":
			compressed_dict[i] = "tar"
		elif i[i.rfind("zip"):] == "zip":
			compressed_dict[i] = "zip"
		else:
			compressed_dict[i] = "no-compression"
	return compressed_dict

def phred_parse (fastqlist, sample_size):
	phred64dict = {}
	counter = 0
	for element in fastqlist:
		switch = False
		for line in open(element):
			if counter > sample_size: break
			if line[0] == "@": continue
			else:
				if switch == False:
					if line[0] == "+":
						switch = True
						counter = counter + 1
					else: continue
				if switch == True:
					if line.find("Z") > -1:
						phred64dict[element] = "64"
						break
				else:
					switch = False
	for element in fastqlist:
		if element not in phred64dict:
			phred64dict[element] = "33"
	return phred64dict

def hypo_dict_parse(fastqlist):
	hypo_dict = {}
	for element in fastqlist:
		if element in hypo_dict:
			continue
		for m in re.finditer('2', element):
			hypothetical = element[:m.start()]+"1"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				break
		for m in re.finditer('R', element):
			hypothetical = element[:m.start()]+"F"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				break
		for m in re.finditer('rev', element):
			hypothetical = element[:m.start()]+"fwd"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				break
		for m in re.finditer('rev', element):
			hypothetical = element[:m.start()]+"fw"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				break
	return hypo_dict

def format_parse(fastq):
	format_dict = {}
	for i in fastq:
		if i[-3:] == "fa" or i[-6:] == ".fasta" :
			format_dict[i] = "fasta"
		else:
			format_dict[i] = "fastq"
	return format_dict

def type_parse(fastq, hypo_dict, mean_read_dict):
	type_dict = {}
	library_size_dict = {}
	for i in fastq:
		if i in type_dict: 
			continue
		if i in hypo_dict:
			type_dict[i] = [1, hypo_dict[i]]
			type_dict[hypo_dict[i]] = [2, i]
		#elif mean_read_dict[i] > 1500:
			# type_dict[i] = ["pb", "no_partner"]			
		else: 
			# type_dict[i] = ["s", "no_partner"]
			type_dict[i] = ["--12", "no_partner"]
	
	for i in fastq:
		if type_dict[i] == 1 or type_dict[i] == 2:
			library_size_dict[i] = (os.stat(i).st_size)*2
		else:
			library_size_dict[i] = (os.stat(i).st_size)
	return type_dict, library_size_dict

# this function creates the prepared_libraries.txt file
def preparation(initial_fastq, sample_size, output_report):
	fastq = remove_false_files(initial_fastq)
	mean_read_dict = {}
	compressed_dict = compression_parse(fastq)
	for i in fastq:
		mean_read_dict[i] = get_mean_read_len(i, sample_size, compressed_dict)

	hypo_dict = hypo_dict_parse(fastq)
	type_dict, library_size_dict = type_parse(fastq, hypo_dict, mean_read_dict)

	phred64dict = phred_parse(fastq, sample_size)
	format_dict = format_parse(fastq)

	report = open(output_report, 'w')
	for i in fastq:
		print (i + "\t" + str(mean_read_dict[i][0]) + "\t" + str(mean_read_dict[i][1]) + "\t" + str(library_size_dict[i]) + "\t" + str(phred64dict[i]) + "\t" + str(type_dict[i][0]) + "\t" + str(type_dict[i][1]) + "\t" + format_dict[i] + "\t" + compressed_dict[i]+"\n")
		report.write(i + "\t" + str(mean_read_dict[i][0]) + "\t" + str(mean_read_dict[i][1]) + "\t" + str(library_size_dict[i]) + "\t" + str(phred64dict[i]) + "\t" + str(type_dict[i][0]) + "\t" + str(type_dict[i][1]) + "\t" + format_dict[i] + "\t" + compressed_dict[i]+"\n")
	report.close()
