# CodonMuSe
**_CodonMuSe_** is a method for analysing the factors that affect codon bias both at the genome-wide and individual gene-level. 

**_CodonMuSe_** estimates
1) The strength of selection acting on nucleotide biosynthesis cost.
2) The strength of selection acting on translational efficiency.
3) The proportion of codon bias attributable to mutational processes.

**_CodonMuSe_** determines the extent to which individual gene sequences are
1) Optimised for cost efficiency.
2) Optimised for translational efficiency.
3) Jointly optimised for both cost and translational efficiency.

## How _CodonMuSe_ works
**_CodonMuSe_** implements the SK model for infering the impact on **Codon** bias of **Mu**tation and **Se**lection acting on nucleotide cost and tranlsational efficiency. The SK model was first described in:

**Seward EA, Kelly S (2016)** Dietary nitrogen alters codon bias and genome composition in parasitic microorganisms. **_Genome Biology_** 17(1):226.

The **_CodonMuSe_** implementation of the SK model was first published in:

**Seward EA, Kelly S (in prep.)** The evolutionary economics of a gene, bacteria balance the cost and efficiency of mRNA.

## Installing _CodonMuSe_
**_CodonMuSe_** is written in python and requires **scipy**. Up-to-date instructions on how to install **scipy** are provided here: http://www.scipy.org/install.html. Once **scipy** is installed, the **_CodonMuSe_** source code can be downloaded from this repositry by clicking the _Clone or download_ link at the top of the page and executed directly as described below.


## Running _CodonMuSe_

`python CodonMuSe.py [OPTIONS] -f <SEQUENCE FILE> -tscan <tRNAscan FILE>`

	-f <FILE>	A FASTA file of protein coding nucleotide seqeunces
	-tscan <FILE>	A tRNA copy number file produced by tRNAscan
	-tc <INT>	The NCBI genetic code identifier goo.gl/ByQOau (Default = 1)
	-ind		Analyse individual genes in adition to a genomewide analysis
	-fix_mb		Fix mutation bias to genome-wide value for individual genes
	-par 		Determine the cost and efficiency optimality of individual genes

## Input Files

Two input files are reqiuerd to run **_CodonMuSe_** 

**1) A FASTA file containing nucleotide sequences for protein coding genes**

Typically this is the full set of protein coding genes from a single species. However, the method can be applied to sequences obtain from transcriptome or metagenome sequencing. By default CodonMuSe only analyses nucleotide sequences that begin with a start codon, end with stop codon, and do not contain in-frame stop codons. Sequences not meeting all of these criteria are ignored and reported in an error file.

**2) tRNA copy number file generated by tRNAscan**

This is a text file containing the output generated by running tRNAscan on the complete genome of the species in question. If you do not include this optional file then CodonMuSe cannot infer the strength of selection acting on codon translational efficiency or determine the extent to which gene sequences have been optimised for translational efficiency.

**Schattner P, Brooks AN, Lowe TM (2005)** The tRNAscan-SE, snoscan and snoGPS web servers for the detection of tRNAs and snoRNAs. **_Nucleic Acids Res_** 33:686–689


## Results files

\<SEQUENCE FILE\>_GenomeWideResults.txt

\<SEQUENCE FILE\>_IndividualGenesResults.txt

\<SEQUENCE FILE\>_OptimisationResults.txt


## Example Dataset

An example dataset is provided. To run CodonMuSe on this dataset execute the following command

`python CodonMuSe.py -f M_pneumoniae_CDS.fasta -tscan M_pneumoniae_tRNAscan.txt -tc 1 -ind -fix_mb -par`
