#!/bin/python
import sys, numpy, os.path, re
import argparse
from Bio import SeqIO

def trimming (library_file, path, commands, job_output, output, remove_originals):
	paired_list = []
	single_list = []
	pacbio_list = []
	backstring = ''
	trimmo_exec = ''

	if len(commands) == 0:
		commands = 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

	for i in os.listdir(path):
		if i[-4:] == '.jar' and i.find('trimmomatic') > -1:
			trimmo_exec = path + i

	if "All_adapters.fa" not in os.listdir(path+"adapters"):
		allada = open(path+"adapters/All_adapters.fa", "w")
		for i in os.listdir(path+"adapters"):
			for line in i:
				allada.write(line[:-1])
			allada.write("\n")
		allada.close()

	for i in open(library_file):
		chunk = i.split()
		if chunk[5] == "1":
			paired_list.append([chunk[0], chunk[6], chunk[4]])
		elif chunk[5] == "2": continue
		elif chunk[5] == "s":
			single_list.append(chunk[0])
		elif chunk[5] == "pb":
			pacbio_list.append(chunk[0])
		else: continue
	
	output_string = ''
	to_remove = []

	for i in paired_list:
		if i[0][0] == '.': i[0] = i[0][1:]
		if i[1][0] == '.': i[1] = i[1][1:]
		output_string = output_string + "java -jar " + trimmo_exec + " PE -phred" + str(i[2]) + " " + i[0] + " " + i[1] + " "+output+"parsed_paired_"+ \
		i[0][i[0].rfind("/")+1:] + " " + output+"parsed_unpaired_"+i[0][i[0].rfind("/")+1:] + " " +output+"parsed_paired_"+i[1][i[1].rfind("/")+1:]\
		 + " " + output+"parsed_unpaired_"+i[1][i[1].rfind("/")+1:] + " " +"ILLUMINACLIP:" + path+"adapters/All_adapters.fa" + ":2:30:10 "+commands+"\n"
		if remove_originals == True:
			to_remove.append(i[0])
			to_remove.append(i[1])

	for i in single_list:
		if i[0] == '.': i = i[1:]
		output_string = output_string + "java -jar " + trimmo_exec + " SE -phred" + str(i[2]) + " " + i[0] + " " + i[1] + " " +output+i[0][:i[0].rfind("/")+1]+"parsed_"+\
		" ILLUMINACLIP:" + path+"adapters/All_adapters.fa" + ":2:30:10 "+commands+"\n"
		if remove_originals == True:
			to_remove.append(i[0])

	for i in to_remove:
		output_string = output_string + "rm " + i + "\n"
	output_file = open(job_output, 'w')
	output_file.write(output_string)
	output_file.close()


