import argparse, os

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--karyon', default='')
	parser.add_argument('--redundans', default='')
	parser.add_argument('--BWA', default='')
	parser.add_argument('--GATK', default='')
	parser.add_argument('--samtools', default='')
	parser.add_argument('--bcftools', default='')
	parser.add_argument('--picardtools', default='')
	parser.add_argument('--SPAdes', default='')
	parser.add_argument('--KAT', default='')
	parser.add_argument('--nQuire', default='')
	parser.add_argument('--SOAPdenovo', default='')
	parser.add_argument('--trimmomatic', default='')
	parser.add_argument('--BUSCO', default='')
<<<<<<< HEAD
	parser.add_argument('--Platanus', default='')
=======
	parser.add_argument('--platanus', default='')
>>>>>>> 51a88c400f6f697c6471911f3d2f1f29b54f90e0
	parser.add_argument('-o', '--output')
	args = parser.parse_args()
	
	
config_file = open(args.output, 'w')

config_file.write('#The config file is parsed by some scripts in the pipeline and used to set basic parameters in the different programs of the pipeline.\n#By default, the program will use this file. The option to select a different configuration file is provided.\n#The file is read by python as simple lines, using special symbols for encoding the necessary information\n#The special symbols and its meaning are the following\n#	With "#" at the start of the line you can add any comment. Any instance of a "#" will be interpreted as a comment, and thus everything written after that will be ignored.\n#	With "+" at the start of the line you indicate the internal name of the program you are going to refer. In cases in which the same program must be called more than once with different options, we suggest independent names for each set of options\n#With "@" at the start of the line you indicate the location of the program you want to call. @ is mandatory, but if empty it will assume that the program is in bashrc\n#	With ">" at the start of the line you select the command line options. The configuration file must be used to store fixed values in the parameters. You can either add everything in a single line, or use several independent lines.\n#With "?" at the start of the line you select the priority of the program. Depending on the mode of the pipeline, some decision making may be necessary. In those cases, the pipeline will use this number to chose, preferring the lowest number.\n#The script parses everything in the form of dictionary. The key are every value encoded with "+". The rest of the information is encoded in a list. Element 0 is a string containing the location of the binary file that must be called. Element 1 is a string containing all the parameters. Element 2 is an integer encoding the priority of the program.\n#If you want to use other parameters we recommend you to copy this file in other location and modify as you please. Naturally, we cannot promise you that the whole pipeline will work if you use a modified version of the configuration file, so try to be sure of any change you make\n\n\n')

config_file.write('+karyon\n')
if args.karyon == "":
	config_file.write('@'+"\n")
else:
	config_file.write('@'+ os.path.abspath(args.karyon)+"/\n")

config_file.write('+redundans\n')
if args.redundans == "":
	config_file.write('@'+"\n")
else:
	config_file.write('@'+ os.path.abspath(args.redundans)+"/\n")

config_file.write('+BWA\n'+"\n")
if args.BWA == "":
	config_file.write('@'+"\n")
else:
	config_file.write('@'+ os.path.abspath(args.BWA)+"/\n")

config_file.write('+GATK\n')
if args.GATK == "":
	config_file.write('@'+"gatk"+"\n")
else:
	config_file.write('@'+ os.path.abspath(args.GATK)+"/gatk"+"\n")
	
config_file.write('+samtools\n')
if args.samtools == "":
	config_file.write('@'+"\n")
else:
	config_file.write('@'+ os.path.abspath(args.samtools)+"/\n")

config_file.write('+bcftools\n')
if args.bcftools == "":
	config_file.write('@'+"\n")
else:
	config_file.write('@'+ os.path.abspath(args.bcftools)+"/\n")
	
config_file.write('+picard-tools\n')
if args.picardtools == "":
	config_file.write('@'+"\n")
else:
	config_file.write('@'+ os.path.abspath(args.picardtools)+"/\n")

config_file.write('+SPAdes\n')
if args.SPAdes == "":
	config_file.write('@')
else:
	config_file.write('@'+ os.path.abspath(args.SPAdes)+"/\n")

config_file.write('+KAT\n')
if args.KAT == "":
	config_file.write('@\n')
else:
	config_file.write('@'+ os.path.abspath(args.KAT)+"\n")

config_file.write('+nQuire\n')
if args.nQuire == "":
	config_file.write('@'+"nQuire\n")
else:
	config_file.write('@'+ os.path.abspath(args.nQuire)+"/nQuire\n")

config_file.write('+SOAPdeNovo\n')
if args.SOAPdenovo == "":
	config_file.write('@\n')
else:
	config_file.write('@'+ os.path.abspath(args.SOAPdenovo)+'/SOAPdenovo-127mer\n')
config_file.write('> -K 63\n')
config_file.write('> -R\n')

config_file.write('+trimmomatic\n')
if args.trimmomatic == "":
	config_file.write('@\n')
else:
	config_file.write('@'+ os.path.abspath(args.trimmomatic)+'/\n')

config_file.write('+BUSCO\n')
if args.BUSCO == "":
	config_file.write('@\n')
else:
	config_file.write('@'+ os.path.abspath(args.BUSCO)+'/\n')
config_file.write('> --auto-lineage-euk\n')
config_file.write('> -m geno\n')
config_file.write('> -f \n')

<<<<<<< HEAD
config_file.write('+Platanus\n')
if args.trimmomatic == "":
	config_file.write('@\n')
else:
	config_file.write('@'+ os.path.abspath(args.Platanus)+'/\n')
=======
config_file.write('+platanus\n')
if args.platanus == "":
	config_file.write('@\n')
else:
	config_file.write('@'+ os.path.abspath(args.platanus)+'/\n')


config_file.close()
