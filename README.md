# CodonMuSel
CodonMuSel is a tool for analysing codon bias both at the genome-wide level and at the gene-level. It estimates 

1) Estimates the strength of selection acting on nucleotide biosynthesis cost.
2) Estimates the strength of selection acting on translational efficiency.
3) Estimates the proportion of codon bias attributable to mutational processes.
4) Determines the extent to which individual gene sequences are optimised for cost efficiency.
5) Determines the extent to which individual gene sequences are optimised for translational efficiency.
6) Determines the extent to which individual gene sequences are jointly optimised for both cost and translational efficiency.

For more details please see:

Seward EA, Kelly S (2016) Dietary nitrogen alters codon bias and genome composition in parasitic microorganisms. Genome Biol 17(1):226.

AND

Seward EA, Kelly S (2016) Dietary nitrogen alters codon bias and genome composition in parasitic microorganisms. Genome Biol 17(1):226.

Input requirements:

Coding sequences for your species of interest in fasta format. (eg. Genus_species.fasta)
Ideally the fasta file should only contain genes that start with a start codon, stop with a stop codon (as codon bias differs along the length of a gene) and are divisable by 3. 
Sequences not meeting these criteria are ignored.

Optional:

A tRNAscan file for your species of interest. (eg. Genus_species_tRNAscan.txt)

(without the tRNAscan file the model will not use translational efficiency as a parameter)

Schattner P, Brooks AN, Lowe TM (2005) The tRNAscan-SE, snoscan and snoGPS web servers for the detection of tRNAs and snoRNAs. Nucleic Acids Res 33(SUPPL. 2):686–689

***Command to run script:***

python Codon_MuSel.py -f Genus_species.fasta -tscan Genus_species_tRNAscan.txt -tc 1 -ind -fix_mb -par

	-f is the fasta file of DNA sequences ***This is not optional***
	-tscan is the species tRNAscan file
	-tc is the translation code you want (standard = 1)
	-ind is if you want to run the script on individual sequences
	-fix_mb is to fix the mutation bias when running on individual species
	-par is to find the pareto optima of individual genes considering cost and translational efficiency

Output:

Species_results_Bestmodelparameters.txt
Species_results_file_individual_genes_Bestmodel_parameters.txt
Species_Paerto_optimisation_results.txt

An example dataset is provided. To run it for this set use the command:

python Codon_MuSel.py -f Mycoplasma_pneumoniae.fasta -tscan Mycoplasma_pneumoniae_tRNAscan.txt -tc 1 -ind -fix_mb -par
