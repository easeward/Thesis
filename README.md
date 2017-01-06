# Codon_bias
This script provides a fast tool for analysing codon bias both at the genome-wide level and at the gene-level.

For more details please see:

Seward EA, Kelly S (2016) Dietary nitrogen alters codon bias and genome composition in parasitic microorganisms. Genome Biol 17(1):226.

Input requirements:

Coding sequences for your species of interest in fasta format. (eg. Genus_species.fasta)

Ideally the fasta file should only contain genes that start with a start codon, stop with a stop codon (as codon bias differs along the length of a gene) and are divisable by 3.

***IT IS IMPORTANT TO ONLY INCLUDE SEQUENCES THAT ARE IN FRAME.***

Sequences not in frame (ie not divisable by 3 with the first nucleotide of the sequence being the first position of the codon) will be wrongly interpreted.

Optional:
A tRNAscan file for your species of interest. (eg. Genus_species_tRNAscan.txt)

(without the tRNAscan file the model will not use translational efficiency as a parameter)
Schattner P, Brooks AN, Lowe TM (2005) The tRNAscan-SE, snoscan and snoGPS web servers for the detection of tRNAs and snoRNAs. Nucleic Acids Res 33(SUPPL. 2):686–689

Command to run script:
python Codon_bias_analysis.py Genus_species.fasta Genus_species_tRNAscan.txt

Output:
Results_Genus_species.txt (This contains the best-model and parameters for your input file).

Additional options:
To run the script on individual genes within a file as well as for the file as a whole add -ind to the command. 
ie. python Codon_bias_analysis.py Genus_species.fasta Genus_species_tRNAscan.txt -ind
It is advisable to fix mutation bias for the organism to the genome-wide level when running each gene individually.
ie. python Codon_bias_analysis.py Genus_species.fasta Genus_species_tRNAscan.txt -ind -fix_mb

Output:
Genus_species_results_file_individual_genes.txt (This contains the best-model and parameters for each sequence in your file).

An example dataset is provided. To run it for this set use the command:
python Mycoplasma_pneumoniae.fasta Mycoplasma_pneumoniae_tRNAscan.txt -ind -fix_mb


