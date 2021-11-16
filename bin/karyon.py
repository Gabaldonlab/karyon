
desc="""Karyon pipeline.
More info at: https://github.com/Gabaldonlab/karyon
"""
epilog="""Author: Miguel A. Naranjo-Ortiz (miguelangelnaranjoortiz@gmail.com) Worcester MA, 04/Nov/2021"""
import sys, os, re
import argparse
import psutil
import pysam
import pandas as pd
import string
import random
from spades_recipee import call_SPAdes
from prepare_libraries import preparation
from trimming_libraries import trimming
from varcall_recipee import var_call
from datetime import datetime

parser = argparse.ArgumentParser(description=desc, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-d', '--output_directory', required=True, help='Directory where all the output files will be generated. Required.')
parser.add_argument('-o', '--output_name', default=False, help='Prefix name for all the output files. If omitted, it will generate a random string. This random string will be the same as the identifier for intermediate files.')
parser.add_argument('-l', '--libraries', required=True, nargs='+', help="Fastq libraries to use for assembly and variant calling. Unsuitable libraries for any of the steps will be ignored. Required.")
parser.add_argument('-F', '--favourite', default=False, help='Sets one library as the prefered one for the variant calling analysis. Otherwise, karyon will select the largest library for performing the variant calling protocol.')
parser.add_argument('-c', '--configuration', default=False, help="Configuration file. By default will use ./configuration.txt as the configuration file.")
parser.add_argument('-g', '--genome_assembler', default="dipspades", choices=['dipspades','dipSPAdes','spades', 'SPAdes','platanus','Platanus', 'soapdenovo', 'SOAPdenovo'], help="Genome assembly software to use. By default it will use dipSPAdes. Options are: dipSPADEs, SPAdes, SOAPdenovo or Platanus.")
parser.add_argument('-T', '--no_trimming', action='store_true', default=False, help='If this tag is active, the program will skip the trimming step.')
parser.add_argument('-A', '--no_assembly', default=False, help='If this tag is active it will skip the assembly step. It requires a reference assembly.')
parser.add_argument('-R', '--no_reduction', action='store_true', default=False, help='If this tag is active, the program will not launch the reduction step of redundans. Remember that the step is used to perform many downstream analyses. If you skip it, the analyses may not make much sense.')
parser.add_argument('-V', '--no_varcall', nargs='+', default=False, help="If this tag is active, the program will skip the variant calling step. Many downstream analyses require this and won't be possible if you skip it.")
parser.add_argument('-B', '--no_busco', default=False, action='store_true', help='If this tag is active, BUSCO analysis will be ommited.')
parser.add_argument('-P', '--no_plot', action='store_true', default=False, help="If this tag is active, the program will omit the plots at the end of the the variant calling step.")
parser.add_argument('-w', '--window_size', default=1000, help='Window size used for some of the analyses. Default is 1000 (1Kb)')
parser.add_argument('-x', '--max_scaf2plot', default=20, help="Maximum number of scaffolds to plot for scaffold-specific plots. Default is 20.")
parser.add_argument('-s', '--scafminsize', default=False, help="Will ignore scaffolds with length below the given threshold")
parser.add_argument('-S', '--scafmaxsize', default=False, help="Will ignore scaffolds with length above the given threshold")
parser.add_argument('-a', '--try_again', default=False, action='store_true', help='Use previous karyon results and skips already computed steps.')
parser.add_argument('-K', '--keep_tmp', action='store_true', default=False, help='If this tag is active, the program will not remove all intermediary files in the folder tmp after it has finished')
parser.add_argument('-i', '--job_id', default=False, help='Identifier of the intermediate files generated by the different programs. If false, the program will assign a name consisting of a string of 6 random alphanumeric characters.')
parser.add_argument('-m', '--memory_limit', default=False, help='Memory limit for all the programs set in Gb. By default it will try to use all memory available.')
parser.add_argument('-M', '--memory_fraction', default=1, help='Proportion of total memory to use by all programs. By default it will use all available memory (default=1), but it may be useful to reduce the percent to avoid freezing other tasks of the computer during peaks.')
parser.add_argument('-n', '--nodes', default=False, help='Number of computation nodes to use. If set a number higher than total, it will use total. If set a number lower than total, it will calculate memory usage based on the fraction of nodes set with respect to total existing nodes.')
args = parser.parse_args()


def id_generator(size=6, chars=string.ascii_uppercase + string.digits): 
	return ''.join(random.choice(chars) for _ in range(size))

###Parses the config file in order to check the parameters of all the programs.###
def parse_config(config):
	config_dict = {}
	prev = 0
	for line in open(config):
		if line[0] == "#": continue
		elif line[0] == "+":
			prev = line[1:-1]
			config_dict[prev] = ["","",""]
		elif line[0] == "@":
			if config_dict[prev][0] != "": continue
			config_dict[prev][0] = line[1:-1]
		elif line[0] == ">":
			config_dict[prev][1] = config_dict[prev][1] + line[1:-1] + " "
		elif line[0] == "?":
			if config_dict[prev][2] != "": continue
			config_dict[prev][2] = line[1:-1] + " "
	return config_dict

###Selects the main library to use. This is set to accelerate the assembly process and improve the results
def select_champion(fastq):
	parse_dict = {}
	for i in open(fastq):
		chunk = i.split()
		if chunk[5] == "2": continue
		else:
			parse_dict[chunk[0]] = chunk[1:]
	champion=[0,'']
	if args.favourite == False:
		for element in parse_dict:
			if int(parse_dict[element][2]) > champion[0]:
				champion = [int(parse_dict[element][2]), element]
	else:
		champion = [0,args.favourite]
	return champion

def exit_program(message):	
	sys.stderr.write("\n%s\n\n"%message)
	sys.exit(1)


def main():
    

	###Defines the location of configuration.txt if setting by default###
	config_path = args.configuration
	if not args.configuration:
		selfpath = os.path.dirname(os.path.realpath(sys.argv[0]))
		config_path = selfpath[:selfpath.rfind('/')]
		config_path = selfpath[:selfpath.rfind('/')]+"/configuration.txt"
	
	true_output = os.path.abspath(args.output_directory)
	if true_output[-1] != "/":
		true_output=true_output+"/"

	###Sets RAM usage options###
	total_nodes = n_nodes = psutil.cpu_count()
	if args.nodes and int(args.nodes) < total_nodes:
		n_nodes = int(args.nodes)

	if not args.memory_limit:
		ram_limit = int(psutil.virtual_memory()[0]/1000000000 * float(args.memory_fraction))
		if n_nodes < total_nodes:
			ram_limit = int(psutil.virtual_memory()[0]/1000000000 * float(args.memory_fraction) * (float(n_nodes)/total_nodes))
	else:
		ram_limit = args.memory_limit * int(args.memory_fraction)
	counter = int(args.max_scaf2plot)
	
	###Sets the job ID and the prefix name for the job. If job ID is not user defined, it produces a random 6 character string. If prefix name is not defined, it uses job ID### 
	job_ID = args.job_id if args.job_id else id_generator()
	name = args.output_name if args.output_name else job_ID

	print ('###############')
	print ('Config. path: '+str(config_path))
	print ("RAM Limit: "+str(ram_limit)+"Gb")
	print ("Nodes: "+str(n_nodes))
	print ("Job ID: "+str(job_ID))
	print ("Job name: "+str(name))
	print ('###############')

	config_dict = parse_config(config_path)
	home = config_dict["karyon"][0]
	if home[-1] != "/": home = home + "/"
	prepared_libs = home + "tmp/" + job_ID + "/prepared_libraries.txt"
	
	path_tmp_jobid = os.path.join(home, "tmp", job_ID)
	if not os.path.exists(os.path.join(home, "tmp")):
		os.mkdir(os.path.join(home, "tmp"))
	prepared_libs = os.path.join(path_tmp_jobid, "prepared_libraries.txt")

	###Checks that the output is not a file. If it does not exist, it creates it.###
	if not os.path.isdir(args.output_directory):
		if os.path.isfile == True: 
			message = "Path is a file" #Should raise an exception an exit the program
			exit_program(message)
		else:
			os.mkdir(args.output_directory)
	elif args.try_again == False:
		os.rmdir(args.output_directory)
		os.mkdir(args.output_directory)
	
	os.mkdir(path_tmp_jobid)


	from karyonplots import katplot, allplots
	katplot("", "", config_dict["KAT"][0], "")


	###Parses the libraries and checks their parameters for downstream analyses. Also performs trimming.###
	print ('###############')
	print ('Preparing libraries')	
	print ('###############')
	libs = ''
	for i in args.libraries:
		libs = libs + " " + os.path.abspath(i)
	preparation(libs.split(), 10000, prepared_libs)

	libs_parsed = ''
	if not args.no_trimming:
		print ('###############')
		print ('Trimmomatic')
		print ('###############')
	
		if config_dict['trimmomatic'][1] == '':
			trimmo_commands = ''
		else:
			trimmo_commands = " -c " + config_dict['trimmomatic'][1]
		trimming(prepared_libs, config_dict["trimmomatic"][0], trimmo_commands, home + "tmp/"+job_ID+"/trimmomatic.job", true_output, False)
		os.system("bash " + home + "tmp/"+job_ID+"/trimmomatic.job")
	
		for i in os.listdir(args.output_directory):
			if i.find("parsed_") > -1:
				libs_parsed = libs_parsed + " " + true_output + i		
		preparation(libs_parsed.split(), 10000, prepared_libs)

	###Parsing library names, including putting absolute paths #
	libstring = ''
	backstring = ''
	for i in open(prepared_libs):
		chunk = i.split()
		if chunk[5] == "1":
			libstring = libstring + os.path.abspath(chunk[0]) + " " + os.path.abspath(chunk[6]) + " "
		elif chunk[5] == "2": continue
		else: backstring = backstring + os.path.abspath(chunk[0]) + " "
	libstring = libstring + backstring

	champion = select_champion(prepared_libs)

	print ('###############')
	print ('Params')
	print ('###############')
	print (args.window_size)
	print (true_output+name+".raw.vcf")
	print (true_output+"redundans_output/scaffolds.filled.fa")
	print (true_output+name+".sorted.bam")
	print (true_output+name+".mpileup")
	print (champion[-1])
	print (config_dict['nQuire'][0])
	print (config_dict["KAT"][0])
	print (home + "tmp/"+job_ID+"/")
	print (true_output)
	print (counter)
	print ('###############')
	
	###Calling spades_recipee.py to generate the assembly job. In the future it should use config file to select the assembly program to use###
	karyonjobfile = open(true_output+name+"_karyon.job", 'a')
	karyonjobfile.write("\n")
	switch = False
	if args.no_assembly == False:
		if args.genome_assembler == "dipspades" or args.genome_assembler == 'dipSPAdes':
			call_SPAdes(prepared_libs, config_dict['SPAdes'][0], true_output, name, config_dict['SPAdes'][1], False, ram_limit, n_nodes)
			assembly = true_output+"dipspades/consensus_contigs.fasta"
		elif args.genome_assembler == "spades" or args.genome_assembler == 'SPAdes':
			call_SPAdes(prepared_libs, config_dict['SPAdes'][0], true_output, name, config_dict['SPAdes'][1], True, ram_limit, n_nodes)
			assembly = true_output+"spades/scaffolds.fasta"
		elif args.genome_assembler == "platanus" or args.genome_assembler == "Platanus":
			if args.no_reduction == True:
				karyonjobfile.write("python2 "+config_dict['redundans'][0]+"redundans.py"+ " -o "+true_output+"redundans_output -i "+libstring+" -t "+str(n_nodes)+" "+config_dict["redundans"][1] + " --noreduction")
			else:
				karyonjobfile.write("python2 "+config_dict['redundans'][0]+"redundans.py"+ " -o "+true_output+"redundans_output -i "+libstring+" -t "+str(n_nodes)+" "+config_dict["redundans"][1])
			assembly = true_output+"redundans_output/scaffolds.filled.fa"
			switch = True
		elif args.genome_assembler == "soapdenovo" or args.genome_assembler == "SOAPdenovo":
			from soap_recipee import soap_recipee
			soap_recipee(prepared_libs, name, true_output+"soapdenovo/", config_dict['SOAPdeNovo'][1], karyonjobfile, config_dict['SOAPdeNovo'][0])
			print ("python3 "+os.path.dirname(__file__)+"/soap_recipee.py -r "+prepared_libs+" -n "+name+" -o "+true_output+"soapdenovo "+ "-j "+true_output+name+"_karyon.job")
			os.system("python3 "+os.path.dirname(__file__)+"/soap_recipee.py -r "+prepared_libs+" -n "+name+" -o "+true_output+"soapdenovo "+ "-j "+true_output+name+"_karyon.job")
			assembly = true_output+"soapdenovo/"+name+".scafSeq"
		else:
			pass
	else:
		assembly = args.no_assembly
	if args.no_reduction == False and switch == False:
		karyonjobfile.write("python2 "+config_dict['redundans'][0]+"redundans.py"+" -f "+ assembly + " -o "+true_output+"redundans_output -i "+libstring+" -t "+str(n_nodes)+" "+config_dict["redundans"][1])
		reduced_assembly = true_output+"redundans_output/scaffolds.filled.fa"
	elif args.no_reduction == False and switch == True:
		reduced_assembly = assembly
	else:
		reduced_assembly = assembly
	busco_options = ""
	if args.no_busco == False:
		for i in config_dict['BUSCO'][1:]:
			busco_options = busco_options + " " + i[:-1]
		karyonjobfile.write("\n")
		karyonjobfile.write(config_dict['BUSCO'][0]+"busco " + "-i " + reduced_assembly + " -o " + name + busco_options + "\n")
		karyonjobfile.write("mv " + name + " " + true_output+name+"_busco\n")
		karyonjobfile.write("mv " + true_output+name+"_busco/short_summary.specific.*.txt " + true_output+name+".busco\n")
		karyonjobfile.write("rm -r busco_downloads\n")
	karyonjobfile.close()

	#5) Create job file that calls all the programs
	if args.no_varcall == False:
		var_call(prepared_libs, config_dict, true_output, name, args.favourite, home, str(ram_limit), str(n_nodes), reduced_assembly)
	os.system ("bash "+true_output+name+"_karyon.job")
	counter = int(args.max_scaf2plot)
	
	def parse_no_varcall(no_varcall):
		vcf, bam, mpileup = '', '', ''
		for i in no_varcall:
			if i[-4:] == ".bam":
				bam = i
			if os.path.isfile(i+".bam") == True:
				bam = i+".bam"
			if os.path.isfile(i+".sorted.bam") == True:
				bam = i+".sorted.bam"
			if i.find("pileup") > -1:
				mpileup = i
			if os.path.isfile(i+".mpileup") == True:
				mpileup = i+".mpileup"
			if os.path.isfile(i+".pileup") == True:
				mpileup = i+".pileup"
			if i[-4:] == ".vcf":
				vcf = i
			if os.path.isfile(i+".vcf") == True:
				vcf = i+".vcf"
			if os.path.isfile(i+"raw.vcf") == True:
				vcf = i+"raw.vcf"
		return vcf, bam, mpileup

	if args.no_varcall == False:
		from karyonplots import katplot, allplots
		from report import report, ploidy_veredict
		katplot(reduced_assembly, champion[1], config_dict["KAT"][0], true_output)
		df = allplots(int(args.window_size), 
			true_output+name+".raw.vcf", 
			reduced_assembly, 
			true_output+name+".sorted.bam", 
			true_output+name+".mpileup", 
			os.path.abspath(champion[-1]), 
			config_dict['nQuire'][0], 
			config_dict["KAT"][0], 
			home + "tmp/"+job_ID+"/", 
			true_output, 
			counter, 
			job_ID, name, args.scafminsize, args.scafmaxsize, args.no_plot)
	else:
		from karyonplots import katplot, allplots
		katplot(reduced_assembly, champion[1], config_dict["KAT"][0], true_output)
		vcf, bam, mpileup = parse_no_varcall(args.no_varcall)
		df = allplots(int(args.window_size), 
			vcf, 
			reduced_assembly, 
			bam, 
			mpileup, 
			os.path.abspath(champion[-1]), 
			config_dict['nQuire'][0], 
			config_dict["KAT"][0], 
			home + "tmp/"+job_ID+"/", 
			true_output, 
			counter, 
			job_ID, name, args.scafminsize, args.scafmaxsize, df, args.no_plot)
	os.path.mkdir(true_output+"/Report/")
	df2 = ploidy_veredict(df, true_output, name)
	report(true_output, name, args.no_reduction)
	df2.to_csv(true_output+"/Report/report"+name+".csv", index=False)

	###We clean the tmp directory###
	if args.keep_tmp == True:
		existence = open(home + "tmp/" + job_ID + '/_keep_existing_', 'w')
		existence.close()

	print ("Now I'm cleaning tmp...")
	if args.keep_tmp == True:
		print ("...but keeping what you told me...")
	for e in os.listdir(home + "tmp/"):
		for i in os.listdir(home + "tmp/"+e):
			if '_keep_existing_' in os.listdir(home + "tmp/"+e): continue
			else:
				os.remove(home + "tmp/"+e+"/"+i)
		if '_keep_existing_' in os.listdir(home + "tmp/"+e): continue
		else: os.rmdir(home + "tmp/"+e)
	if args.keep_tmp == True:
		print ("... tmp files havee been kept")
	else:	
		print ("... removed tmp files!")

if __name__ == '__main__':
	t0  = datetime.now()
	try:
		main()
	except KeyboardInterrupt:
		sys.stderr.write("\n Ctrl-C pressed!      \n")
	dt  = datetime.now()-t0
	sys.stderr.write("#Time elapsed: %s\n" % dt)

