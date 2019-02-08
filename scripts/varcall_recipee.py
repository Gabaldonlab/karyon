#!/bin/python
import sys, numpy, os.path, re
import argparse
from Bio import SeqIO
'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--libraries', required=True, nargs='+', help="Fastq libraries to use for assembly and variant calling. Unsuitable libraries for any of the steps will be ignored")
    parser.add_argument('-n', '--name',  required=True, help="Prefix to use with the output files")
    parser.add_argument('-o', '--output', required=True, help="Output directory")

args = parser.parse_args()

if args.output[-1] == "/":
	args.output = args.output[:-1]
if os.path.isdir(args.output) == False:
	os.mkdir(args.output)

locspp = args.output + "/" + args.name
fastq = args.libraries


libraries_location=args.libraries[0][:fastq[0].rfind("/")]
location=''
'''
def select_champion(fastq, favourite):
	parse_dict = {}
	for i in open(fastq):
		chunk = i.split()
		if chunk[5] == "2": continue
		else:
			parse_dict[chunk[0]] = chunk[1:]
	champion=[0,'']		
	if favourite == False:
		for element in parse_dict:
			if int(parse_dict[element][2]) > champion[0]:
				champion = [int(parse_dict[element][2]), element]
	else:
		champion = [0,args.favourite]
	return champion, parse_dict

def job_description (fastqlist):
	phred64dict = {}
	for element in fastqlist:
		switch = False
		for line in open(element):
			if line[0] == "@": continue
			else:
				if switch == False:
					if line[0] == "+":
						switch = True
					else: continue
				if switch == True:
					if line.find("0") > -1:
						phred64dict[element] = False
						break
				else:
					switch = False
	for element in fastqlist:
		if element not in phred64dict:
			phred64dict[element] = True
	for library in fastq:
		if phred64dict[library] == True:
			SeqIO.convert(library, 'fastq-illumina', library+".converted.fq", 'fastq-sanger')
	return phred64dict


def create_hypo_dict(fastq):
	hypo_dict = {}
	for element in fastq:
		for m in re.finditer('2', element):
			hypothetical = element[:m.start()]+"1"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				fastq.remove(element)
				break
		for m in re.finditer('R', element):
			hypothetical = element[:m.start()]+"F"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				fastq.remove(element)
				break
		for m in re.finditer('rev', element):
			hypothetical = element[:m.start()]+"fwd"+element[m.end():]
			if os.path.exists(hypothetical):
				hypo_dict[hypothetical] = element
				fastq.remove(element)
				break
	return hypo_dict			

#if os.path.exists(args.output+"/"+args.name+".karyon.txt"): os.remove(args.output+"/"+args.name+".karyon.txt")

def var_call(fastq, config_dict, output, name, favourite, home, memory, nodes):
	outputfile = output+name+"_karyon.job"
	parse_dict = {}
	libstring = ' '
	backstring = ''
	for i in open(fastq):
		chunk = i.split()
		if chunk[5] == "2": continue
		parse_dict[chunk[0]] = chunk[1:]
		if chunk[5] == "1":
			libstring = libstring + os.path.abspath(chunk[0]) + " " + os.path.abspath(chunk[6]) + " "	
		elif chunk[5] == 's' : continue
		else: backstring = backstring + os.path.abspath(chunk[5]) + " " + os.path.abspath(chunk[0]) + " "
	libstring = libstring + backstring
	
	bash_job = open(outputfile, "a")
	if favourite == False:
		locspp = output + name
		locspp2 = locspp
	else:
		locspp = output + name
		locspp = output + name + "_" + favourite[:favourite.rfind(".")+1]
	pairs = ''
	bash_job.write("\n")
	bash_job.write("cp "+output+"redundans_output/scaffolds.filled.fa "+locspp+".fasta\n\n")
	bash_job.write("bwa index "+locspp+".fasta\n\n")
	
	champion, parse_dict = select_champion(fastq, favourite)

	bash_job.write("java -Xmx"+memory+"g -jar " + config_dict["picard-tools"][0]+"CreateSequenceDictionary.jar R="+locspp+'.fasta O='+locspp+'.dict\n\n')
	bash_job.write(config_dict["BWA"][0]+"bwa index "+locspp+'.fasta\n\n')
	if parse_dict[champion[1]][4] == '1':
		bash_job.write("python " + home + "scripts/launch_bwa.py -r "+locspp+".fasta -f1 "+os.path.abspath(champion[1])+" -f2 "+os.path.abspath(parse_dict[champion[1]][5])+" -n "+locspp+"\n\n")	
	if parse_dict[champion[1]][4] == 's':
		bash_job.write("python " + home + "scripts/launch_bwa.py -r "+locspp+".fasta -f1 "+os.path.abspath(champion[1])+" -n "+locspp+"\n\n")
	#print "python " + home + "scripts/launch_bwa.py -r "+locspp+".fasta -f1 "+os.path.abspath(champion[1])+" -f2 "+os.path.abspath(parse_dict[champion[1]][5])+" -n "+locspp+"\n\n"
	bash_job.write(config_dict["samtools"][0]+"samtools index "+locspp2+'.sorted.bam\nsamtools faidx '+locspp+'.fasta\n\n')
	bash_job.write(config_dict["samtools"][0]+"samtools index "+locspp+'.sorted.bam\n')
	bash_job.write(config_dict["samtools"][0]+"samtools faidx "+locspp+'.fasta\n')
	bash_job.write(config_dict["GATK"][0] + " --java-options -Xmx"+memory+"G HaplotypeCaller -R "+locspp+'.fasta -I '+locspp+'.sorted.bam -O '+locspp+'.raw.vcf\n\n')
	# bash_job.write(config_dict["samtools"][0]+"samtools mpileup -Q 0 -M 0 -e 0 "+locspp+'.sorted.bam > '+locspp+'.mpileup\n')
	bash_job.write(config_dict["bcftools"][0]+"bcftools mpileup --fasta-ref " + locspp + ".fasta " + config_dict["bcftools"][1]+" "+locspp+'.sorted.bam > '+locspp+'.mpileup\n')	
