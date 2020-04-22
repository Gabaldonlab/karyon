import argparse

if __name__ == '__main__':	
	parser.add_argument('--karyon', default='')
	parser.add_argument('--redundans', default='')
	parser.add_argument('--BWA', default='')
	parser.add_argument('--GATK', default='')
	parser.add_argument('--samtools', default='')
	parser.add_argument('--bcftools', default='')
	parser.add_argument('--picard-tools', default='')
	parser.add_argument('--SPADes', default='')
	parser.add_argument('--KAT', default='')
	parser.add_argument('--nQuire', default='')
	parser.add_argument('--SOAPdeNovo', default='')
	parser.add_argument('-o', '--output')
	args = parser.parse_args()
	
	
config_file = open(args.output, 'w')

config_file.write(#The config file is parsed by some scripts in the pipeline and used to set basic parameters in the different programs of the pipeline.
#By default, the program will use this file. The option to select a different configuration file is provided.
#The file is read by python as simple lines, using special symbols for encoding the necessary information
#The special symbols and its meaning are the following
#	With "#" at the start of the line you can add any comment. Any instance of a "#" will be interpreted as a comment, and thus everything written after that will be ignored.
#	With "+" at the start of the line you indicate the internal name of the program you are going to refer. In cases in which the same program must be called more than once with different options, we suggest independent names for each set of options
#	With "@" at the start of the line you indicate the location of the program you want to call. @ is mandatory, but if empty it will assume that the program is in bashrc
#	With ">" at the start of the line you select the command line options. The configuration file must be used to store fixed values in the parameters. You can either add everything in a single line, or use several independent lines.
#	With "?" at the start of the line you select the priority of the program. Depending on the mode of the pipeline, some decision making may be necessary. In those cases, the pipeline will use this number to chose, preferring the lowest number.

#The script parses everything in the form of dictionary. The key are every value encoded with "+". The rest of the information is encoded in a list. Element 0 is a string containing the location of the binary file that must be called. Element 1 is a string containing all the parameters. Element 2 is an integer encoding the priority of the program.
#If you want to use other parameters we recommend you to copy this file in other location and modify as you please. Naturally, we cannot promise you that the whole thing will work if you use a modified version of the configuration file, so try to be sure of any change you make


)
config_file.write('+karyon\n')
config_file.write('@'+ args.karyon)

config_file.write('+redundans\n')
config_file.write('@'+ args.redundans)

config_file.write('+BWA\n')
config_file.write('@'+ args.BWA)

config_file.write('+GATK\n')
config_file.write('@'+ args.GATK)

config_file.write('+samtools\n')
config_file.write('@'+ args.samtools)

config_file.write('+bcftools\n')
config_file.write('@'+ args.bcftools)

config_file.write('+picard-tools\n')
config_file.write('@'+ args.picard-tools)

config_file.write('+SPADes\n')
config_file.write('@'+ args.SPADes)

config_file.write('+KAT\n')
config_file.write('@'+ args.KAT)

config_file.write('+nQuire\n')
config_file.write('@'+ args.nQuire)

config_file.write('+SOAPdeNovo\n')
config_file.write('@'+ args.SOAPdeNovo)

config_file.close()
