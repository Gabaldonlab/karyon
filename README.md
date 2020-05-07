This repository contains the Karyon pipeline.

Karyon is a pipeline for the assembly and analysis of highly heterozygous genomes. It uses redundans (Pryszcz & Gabaldón, 2016) to reduce heterozygosity during the assembly process, and then maps the original libraries against the reduced assembly to analyze the distribution of heterozygous regions. With this information, it generates a series of plots that can aid researchers to generate informed hypotheses with regard of the architecture of their genomes.

Scripts contained in this repository:
1) karyon.py -> complete pipeline, including genome assembly, assembly reduction, SNP calling and plot generation
2) prepare_libraries.py -> karyon dependency. It uses Trimmomatic () to trim input libraries before genome assembly.
3) spades_recipee.py -> Karyon dependency. It generates a file that launches dipSPAdes () with the input.
4) varcall_recipee.py -> Karyon dependency. It generates a file that launches all steps in the SNP calling pipeline.
5) karyonplots.py -> Karyon dependency. It generates all the plots as part of the Karyon pipeline.
6) all_plots.py -> Standalone version of karyonplots.py. It allows the user to input karyon results to generate the plots again.

## Authors 
* **Miguel Ángel Naranjo Ortiz** - *Pipeline work* - [MANaranjo](https://github.com/MANaranjo)
* **Manuel Molina Marín** - *Docker work* - [manumolina](https://github.com/manumolina)

## License 
This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.