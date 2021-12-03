from FastaIndex import FastaIndex
import sys
import pandas as pd
import numpy as np
import scipy
from Bio import SeqIO
import os

def fasta_stats(fasta):
	faidx = FastaIndex(fname)
	outlist.append(faidx.stats())
	return(outlist)

def get_busco(busco_file):
	my_busco = open(busco_file, 'w')
	switch = False
	buscolist = []
	for line in my_busco:
		if switch == True and len(line) > 0:
			a = line.split(":").split("%")
			buscolist.append(float(a[1]))
			buscolist.append(float(a[3]))
			buscolist.append(float(a[5]))
			buscolist.append(float(a[7]))
			buscolist.append(float(a[9]))
			switch == False
		if line.find("*****") > -1:
			switch = True
		else:
			continue
	return (buscolist)

def get_nQuire(nQuirefile):
	nQlist = []
	for line in open(nQuirefile):
		if line.split()[0] == "file" or len(line) == 0: continue
		else:
			for i in line[:-1].split()[1:5]:
				nQlist.append(float(i))
	return(nQlist)

def get_flagstat(flagstat):
	mapped = 0
	prop_paired = 0
	for line in open(flagstat):
		if line.find("%") > -1 and line.find("mapped") > -1:
			mapped = float(line.split("(")[-1].split("%")[0])
		if line.find("%") > -1 and line.find("properly") > -1:
			prop_paired = float(line.split("(")[-1].split("%")[0])
	return mapped, prop_paired
	
def report(true_output, name, df, no_reduction, no_red_assembly, no_busco, window_size):
	true_output = os.path.abspath(true_output)
	if true_output[-1] != "/":
		true_output=true_output+"/"
	report = open(true_output+"/Report/"+name+"report.txt", "w")
	#fastats = fasta_stats([true_output+name+".fasta"])
	fastats = FastaIndex(true_output+name+".fasta").stats()
	nQlist = get_nQuire(true_output+name+".lrdtest")
	mapped, prop_paired = get_flagstat(true_output+name+".flagstat")
	report.write("###GLOBAL STATS###\n")
	report.write("Scaffolds:\t"+str(fastats[1])+"\n")
	report.write("Assembly size:\t"+str(fastats[2])+"\n")
	report.write("GC%:\t"+str(fastats[3])+"\n")
	if fastats[3] < 35:
		report.write("WARNING: GC content is very low!\n")
	if fastats[3] > 65:
		report.write("WARNING: GC content is very high!\n")
	report.write("Scaffolds > 1000bp:\t"+str(fastats[4])+"\n")
	report.write("Length of scaffolds > 1000bp:\t"+str(fastats[5])+"\n")
	report.write("N50:\t"+str(fastats[6])+"\n")
	report.write("N90:\t"+str(fastats[7])+"\n")
	if no_busco != True:
		buscoval = get_busco(true_output+name+".busco")
		if no_reduction != True:
			buscored = get_busco(true_output+name+"_no_reduc.busco")
		report.write("BUSCO scores:")
		report.write("Complete: "+str(buscoval[0])+"%, of which "+str(buscoval[1])+"% are single copy and "+str(buscoval[2])+"% are duplicated.")
		if buscoval[0] < 90 and buscoval > 70:
			report.write("WARNING: Your genome has a low BUSCO completeness.\n")
		if buscoval <= 70:
			report.write("WARNING: Your genome has a VERY low BUSCO completeness!\n")
		report.write("Fragmented: "+buscoval[3]+"%. Missing: "+buscoval[4])
		if no_red_assembly != False and no_reduction != False:
			no_red_busco = get_busco(true_output+name+"_no_reduc.busco")
			report.write("Values before reduction: C: "+str(no_red_busco[0])+", SC: "+str(no_red_busco[1])+", D: "+str(+no_red_busco[2])+", F:"+str(no_red_busco[3])+", M: "+str(no_red_busco[4])+"\n")
			if buscoval[0]/no_red_busco[0] < 0.9:
				report.write("WARNING: Reduction has caused loss of information, based on BUSCO values before and after reduction.\n")
	#report.write("BUSCO (prokaryotic) scores:")
	report.write("\n")
	report.write("###ALIGNMENT STATS###\n")
	report.write("% of aligned reads:" + str(mapped) + "%\n")
	if mapped < 90:
		report.write("WARNING: The percent of mapped reads is low\n")
	report.write("% of properly paired reads:" + str(prop_paired) + "%\n")
	if mapped < 80:
		report.write("WARNING: The percent of properly paired reads is low\n")
	if no_reduction != True:
		redustats = FastaIndex(no_red_assembly).stats()
		report.write("###REDUCTION STATS###\n")
		report.write("Scaffolds: From "+str(redustats[1])+" to "+str(fastats[1])+" (" +str(100*(redustats[1]-fastastats[1])/redustats[1])+"%)\n")
		report.write("Assembly size: From "+str(redustats[2])+" to "+str(fastats[2])++" (" +str(100*(redustats[1]-fastastats[1])/redustats[1])+"%)\n")
		report.write("Scaffolds > 1000bp:: From "+str(redustats[4])+" to "+str(fastats[4])+" (" +str(100*(redustats[1]-fastastats[1])/redustats[1])+"%)\n")
		report.write("Length of scaffolds > 1000bp:: From "+str(redustats[5])+" to "+str(fastats[5])+" (" +str(100*(redustats[1]-fastastats[1])/redustats[1])+"%)\n")
		report.write("N50: From "+str(redustats[6])+" to "+str(fastats[6])+" (" +str(100*(redustats[1]-fastastats[1])/redustats[1])+"%)\n")
		report.write("N90: From "+str(redustats[7])+" to "+str(fastats[7])+" (" +str(100*(redustats[1]-fastastats[1])/redustats[1])+"%)\n")
	report.write("\n\n")
	report.write("###PLOIDY ANALYSIS###\n")
	report.write("Global nQuire free model score:\t"+str(nQlist[0])+"\n")
	report.write("Global nQuire diploid score:\t"+str(nQlist[1])+"\n")
	report.write("Global nQuire triploid score:\t"+str(nQlist[2])+"\n")
	report.write("Global nQuire tetraploid score:\t"+str(nQlist[3])+"\n")
	ploid_dict = {0:"Haploid, most likely", 1:"Haploid", 2:"Diploid", 3:"Triploid", 4:"Tetraploid"}
	a = df.count()
	b = []
	pdict = df['ploidy'].value_counts().to_dict()
	for i in range(0,5):
		b.append(str(pdict[i]/a))
	report.write("Percent of haploid windows:\t"+b[1]+"%\n")
	report.write("Percent of diploid windows:\t"+b[2]+"%\n")
	report.write("Percent of triploid windows:\t"+b[3]+"%\n")
	report.write("Percent of tetraploid windows:\t"+b[4]+"%\n")
	report.write("Percent of unassigned windows:\t"+b[0]+"%\n")
	report.write("Main ploidy is:\t"+ploid_dict[b.index(max(b))]+"\n")
	report.write("Diploid regions have an average of "+str(df.loc[df["ploidy"] == 2].SNPs.mean()/window_size)+"% of variant loci\n")
	report.write("\n")
	report.write("###Per scaffold analysis###\n")
	fastadict = SeqIO.index(true_output+name+".fasta", "fasta")
	for entry in fastadict:
		contigdata = df["contig" == entry.id]
		nwindows = contigdata.count()
		for i in range(1,4):
			if contigdata["ploidy" == i]/nwindows > 0.5 and i != b.index(max(b)):
				report.write("Scaffold " + fastadict[entry].id + ", with a length of " + str(len(entry.seq))+ "base pairs, has a " + str(contigdata["ploidy" == i]/nwindows) + "% of windows that are" + ploid_dict[i] + ".\n")

def ploidy_veredict(df, true_output, name, window_size):
	b = df.loc[df.diplo_score != np.NAN] 
	a = b.loc[b.diplo_score > b.triplo_score].loc[b.diplo_score > b.triplo_score].loc[:,"mean_cov"]
	empir_mean, empir_stdev = np.mean(a), np.std(a)
	hap_normal, tetra_normal = np.random.normal(empir_mean/2, empir_stdev/2, 10000), np.random.normal(empir_mean*2, empir_stdev*2, 10000)
	hap25, hap75 = np.quantile(hap_normal, 0.25), np.quantile(hap_normal, 0.75)
	tetra25, tetra75 = np.quantile(tetra_normal, 0.25), np.quantile(tetra_normal, 0.75)
	for i in range(0, len(df)):
		if df["mean_cov"][i] > hap25 and df["mean_cov"][i] < hap75 and (df["SNPs"][i])/window_size <= 0.005:
			df.at[i, "ploidy"] = 1
		elif df["diplo_score"][i] != np.NAN:
			if df["mean_cov"][i] > tetra25 and df["mean_cov"][i] < tetra75 and df["tetra_score"][i] > df["diplo_score"][i] and df["tetra_score"][i] > df["triplo_score"][i]:
				df.at[i, "ploidy"] = 4
			elif df["triplo_score"][i] > df["diplo_score"][i]:
				df.at[i, "ploidy"] = 3
			elif df["diplo_score"][i] > df["triplo_score"][i]:
				df.at[i, "ploidy"] = 2
			else:
				df.at[i, "ploidy"] = 0
		else:
			df.at[i, "ploidy"] = 0
	return(df)
		
