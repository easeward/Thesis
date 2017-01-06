# -*- coding: utf-8 -*-

import glob
import sys
import time
import math
from scipy import optimize
import numpy as np
import random
import string

complementary_nuc = { 'A' : 'T', 'T' : 'A', 'C' :'G', 'G' : 'C'}
alphabetic_codons = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
proteins = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 'B': ['TAA', 'TGA', 'TAG'], 'C' : ['TGT', 'TGC'], 'D' : [ 'GAT', 'GAC'], 'E' : ['GAA','GAG'], 'F' : ['TTT', 'TTC'], 'G' : ['GGT', 'GGC', 'GGA', 'GGG'], 'H' : ['CAT', 'CAC'], 'I' : ['ATT', 'ATC', 'ATA'], 'K' : ['AAG', 'AAA'], 'L' : ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'M' : ['ATG'], 'N' : ['AAT', 'AAC'], 'P' : ['CCT', 'CCC', 'CCA', 'CCG'], 'Q' : ['CAA', 'CAG'], 'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],  'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'T' : ['ACT', 'ACC', 'ACA', 'ACG'],  'V' : ['GTT', 'GTC', 'GTA', 'GTG'], 'W' : ['TGG'], 'Y' : ['TAT', 'TAC']}
codon_trans_standard = {'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C', 'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C', 'TTA' : 'L', 'TCA' : 'S', 'TAA' : 'B', 'TGA' : 'B', 'TTG' : 'L', 'TCG' : 'S', 'TAG' : 'B', 'TGG' : 'W', 'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R', 'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R', 'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R', 'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R', 'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S', 'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S', 'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R', 'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R', 'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G', 'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G', 'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G', 'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'}
sij_values = {'AT' : 0, 'CG' : 0, 'GC' : 0, 'TA' : 0, 'TG' : 0.68, 'GT' : 0.561, 'AC' : 0.28, 'CA' : 0.89, 'AA' : 0.9999, 'AG' : 0.9999, 'CC' : 0.9999, 'CT' : 0.9999, 'GG' : 0.9999, 'GA' : 0.9999, 'TT' : 0.9999, 'TC' : 0.9999}
GC_content = {'A' : 0, 'T' : 0, 'C' : 1, 'G' : 1}
N_content = {'A' : 5, 'T' : 2, 'C' : 3, 'G' : 5} 
bases = ['A', 'T', 'C', 'G']

def import_sequence(file):
	species = str(file[:len(file)-6])
	global codon_count
	codon_count = {}
	protein_count = {}
	for key in codon_trans_standard:
		codon_count[key] = 0
	for key in proteins:
		protein_count[key] = 0;
	codon = ""
	GC_percent = 0
	num_codons = 0
	with open(file) as CDSfile:
		for line in CDSfile:
			line = str(line.rstrip().strip())
			if line[0] == ">":
				acc = line[len(species)+2:]
			else:
				p = 1
				for letter in line:
					if p % 3 == 0:
						p += 1
						codon = codon + letter
						if codon in codon_trans_standard:
							codon_count[codon] = codon_count[codon] + 1
							protein_count[codon_trans_standard[codon]] = protein_count[codon_trans_standard[codon]] + 1 #print "%d %s %d %s" %(codon_count[codon], codon, protein_count[codon_trans_standard[codon]], codon_trans_standard[codon])
							GC_percent = GC_percent + GC_content[codon[0]] + GC_content[codon[1]] + GC_content[codon[2]]
							num_codons = num_codons + 1
						elif "N" in codon:
							print "There is an N in sequence %s, that codon is being ignored "%(codon)
						else:
							print "%s is not a standard codon and won't be counted" %(codon)
						codon = ""
					else:
						p += 1
						codon = codon + letter
	prob_GC = "%.4f" %(float(GC_percent)/(num_codons*3))
	global relative_codon_use
	relative_codon_use = {}
	for codon in codon_trans_standard:
		if codon_count[codon] != 0:
			relative_codon_use[codon] = "%.4f" %(float(codon_count[codon])/protein_count[codon_trans_standard[codon]])
		elif protein_count[codon_trans_standard[codon]] !=0:
			relative_codon_use[codon] = 0
			print "The codon %s isn't used, is this odd?\n" %(codon)
	log_likelihood = float(0)
	for codon in relative_codon_use:
		probability = float(relative_codon_use[codon])
		if probability == 0:
			probability = 0.000000000001
		probability = math.log(probability)
		log_likelihood = log_likelihood + probability*codon_count[codon]
	final_log_likelihood = "%.2f" %(log_likelihood)
	#Mutation_bias = Mb	Selection on nucleotide cost = Ns		Selection on Translational Efficiency = Te
	if len(sys.argv) > 2 and "tRNAscan" in sys.argv[2]:
		get_tAI_values()
		Options = ['Mb', 'Ns', 'Te'] # 
	else:
		Options = ['Mb', 'Ns']
	#AIC = 2k - 2ln(L) # minimum AIC is best
	best_log_likelihood = 10000000000000000000
	AIC = 10000000000000000000
	best = ''
	i = 0
	while i < len(Options):
		opt = Options[i]
		model = opt
		number_para = 1
		res_model, log_likelihood_model, r2_model = get_math_function(model, 0, 0, 0)
		AIC_model = (2*number_para)-(2*(log_likelihood_model))
		if AIC_model < AIC:
			best_res, best_log_likelihood, best_r2, best, AIC = res_model, log_likelihood_model, r2_model, model, AIC_model
		j = i + 1
		while j < len(Options):
			opt_2 = Options[j]
			model = opt+"_"+opt_2
			number_para = 2
			res_model, log_likelihood_model, r2_model = get_math_function(model, 0, 0, 0)
			AIC_model = (2*number_para)-(2*(log_likelihood_model))
			if AIC_model < AIC:
				best_res, best_log_likelihood, best_r2, best, AIC = res_model, log_likelihood_model, r2_model, model, AIC_model
			k = j + 1
			while k < len(Options):
				opt_3 = Options[k]
				model = opt+"_"+opt_2+"_"+opt_3
				number_para = 3
				res_model, log_likelihood_model, r2_model = get_math_function(model, 0, 0, 0)
				AIC_model = (2*number_para)-(2*(log_likelihood_model))
				if AIC_model < AIC:
					best_res, best_log_likelihood, best_r2, best, AIC = res_model, log_likelihood_model, r2_model, model, AIC_model
				k = k + 1
			j = j + 1
		i = i + 1
	f1=open("Results_"+species+".txt", "a") #f1=open("Results_"+best+"_for_"+species+".txt", "a")
	out = "Species	Log_likelihood	Best_log_likelihood	AIC	R2	Mutation_bias	Selection_on_Nucleotide_cost	Selection_on_Translational_efficiency"
	f1.write(out)
	if ('Mb' in best):
		Mutation_bias = best_res[0]
	else:
		Mutation_bias = 0
	if ('Ns' in best):
		Nucleotide_cost = best_res[1]
	else:
		Nucleotide_cost = 0
	if ('Te' in best):
		Translational_efficiency = best_res[2]
	else:
		Translational_efficiency = 0
	out = "\n"+species+"\t"+final_log_likelihood+"\t"+str(best_log_likelihood)+"\t"+str(AIC)+"\t"+str(best_r2)+"\t"+str(Mutation_bias)+"\t"+str(Nucleotide_cost)+"\t"+str(Translational_efficiency)
	f1.write(out)
	f1.close()
	return best, best_res[0]
	
def per_gene_analysis(file, best_model):
	species = str(file[:len(file)-6])
	results_name = species+"_results_file_individual_genes.txt"
	f1=open(results_name, "a")
	header = "Accession	Log_likelihood	R2	Mutation_bias	Selection_on_Nucleotide_cost	Selection_on_Translational_efficiency"
	f1.write(header)
	with open(file) as CDSfile:
		if len(sys.argv) > 2 and "tRNAscan" in sys.argv[2]:
			get_tAI_values()	
		for line in CDSfile:
			line = str(line.rstrip().strip())
			if line[0] == ">":
				acc = line[1:]
				f1.write('\n'+acc+'\t')
				global codon_count
				codon_count = {}
				protein_count = {}
				for key in codon_trans_standard:
					codon_count[key] = 0
				for key in proteins:
					protein_count[key] = 0;
				codon = ""
			else:
				p = 1
				for letter in line:
					if p % 3 == 0:
						p += 1
						codon = codon + letter	
						if codon in codon_trans_standard:
							codon_count[codon] = codon_count[codon] + 1
							protein_count[codon_trans_standard[codon]] = protein_count[codon_trans_standard[codon]] + 1
						elif "N" in codon:
							print "There is an N in sequence %s, that codon is being ignored "%(codon)
						else:
							print "%s is not a standard codon and won't be counted" %(codon)
						codon = ""
					else:
						p += 1
						codon = codon + letter
				global relative_codon_use
				relative_codon_use = {}
				for codon in codon_trans_standard:
					if codon_count[codon] != 0:
						relative_codon_use[codon] = "%.4f" %(float(codon_count[codon])/protein_count[codon_trans_standard[codon]])
					elif protein_count[codon_trans_standard[codon]] !=0:
						relative_codon_use[codon] = 0
				likelihood = float(0)
				for codon in relative_codon_use:
					probability = float(relative_codon_use[codon])
					if probability == 0:
						probability = 0.000000000001
					probability = math.log(probability)
					likelihood = likelihood + probability*codon_count[codon]
				final_likelihood = "%.2f" %(likelihood)
				res_model, likelihood_model, r2_model = get_math_function(best_model, 0, 0, 0)
				if fixed_Mb == 'moveable':
					print "Not fixed Mb\n"
				else: 
					res_model[0] = fixed_Mb
				if ('Mb' in best_model):
					Mutation_bias = res_model[0]
				else:
					Mutation_bias = 0
				if ('Ns' in best_model):
					Nucleotide_cost = res_model[1]
				else:
					Nucleotide_cost = 0
				if ('Te' in best_model):
					Translational_efficiency = res_model[2]
				else:
					Translational_efficiency = 0
				out = str(likelihood_model)+"\t"+str(r2_model)+"\t"+str(Mutation_bias)+"\t"+str(Nucleotide_cost)+"\t"+str(Translational_efficiency)
				f1.write(out)
	f1.close()

def get_math_function(model_type, start_Mb, start_Ns, start_Te):
	equati = write_equation(model_type)
	def equation_real(x):
		Mb, Ns, Te = x
		try:
			return eval(equati)
		except OverflowError:
			return 10000000000000000000
		except ZeroDivisionError:
			return 10000000000000000000
		except ValueError:
			return 10000000000000000000
	x0 = np.asarray((start_Mb, start_Ns, start_Te))
	result_parameters = optimize.fmin(equation_real, x0)
	result_log_likelihood = -(equation_real(result_parameters))
	r2_value = get_r2_value(result_parameters, model_type)
	return result_parameters, result_log_likelihood, r2_value

def get_r2_value(x, model_type):
	Mb = x[0]
	Ns= x[1]
	Te = x[2] 
	top_line = 0
	bottom_line = 0
	mean = 0
	number_codons_used = 0
	for codon in relative_codon_use:
		probability_codon = get_probability_codon(codon, model_type, Mb, Ns, Te)
		mean = mean + probability_codon
		number_codons_used = number_codons_used + 1
	mean = mean/number_codons_used
	for codon in relative_codon_use:
		probability_codon = get_probability_codon(codon, model_type, Mb, Ns, Te)
		relative_codon = float(relative_codon_use[codon])
		value = probability_codon - relative_codon
		value = value*value
		top_line = top_line + value
		second_value = probability_codon - mean
		second_value = second_value*second_value
		bottom_line = bottom_line + second_value
	top_over_bottom = float(top_line)/bottom_line
	return (1-top_over_bottom)
	
def get_probability_codon(input_codon, model_type, Mb, Ns, Te):
	prot = codon_trans_standard[input_codon]
	largest_value = 0
	for syn in proteins[prot]:
		if  "Te" in model_type:
			if tAI_value[syn] > largest_value:
				largest_value = tAI_value[syn]
	upper = {}
	for syn in proteins[prot]:
		upper[syn] = ""
	parameters = model_type.split('_')
	number_parameters = len(parameters)
	for parameter in parameters:
		number_parameters -= 1
		for syn in proteins[prot]:
			if "Mb_fixed" in model_type:
				additional = get_upper_equation(syn, parameter, Mb, Ns, Te, 'for_probability', fixed_Mb, 'off')
			else:
				additional = get_upper_equation(syn, parameter, Mb, Ns, Te, 'for_probability', 'not', 'off')
			upper[syn] = upper[syn]+additional
			if number_parameters == 0:
				upper[syn] = upper[syn][:len(upper[syn])-1]
	base = "("
	for syn in proteins[prot]:
		base = base+upper[syn]+" + "
	base = base[:len(base)-3]+")"
	if model_type == "Te" and prot == 'B':
		return float(relative_codon_use[input_codon])
	else:
		to_be_logged = eval(upper[input_codon]+"/"+base)
		if to_be_logged <=0.00000001:
			return 0.00000001
		else:
			return to_be_logged

def get_upper_equation(codon, parameter, Mb, Ns, Te, equation_type, fixed_Mb, shuffle_type):
	if equation_type == "for_probability":
		Mb = '%.6f'%Mb
		Ns = '%.6f'%Ns
		Te = '%.6f'%Te
	if parameter == 'Mb':
		if shuffle_type == 'off':
			GC_codon = '%.4f'%(GC_unshuffled[codon])
		else:
			GC_codon = '%.4f'%(GC_shuffled[codon])
		if fixed_Mb == 'not':
			return "(math.exp("+Mb+"*"+GC_codon+"))*"
		else:
			fixed_Mb = '%.8f'%fixed_Mb
			return "(math.exp("+fixed_Mb+"*"+GC_codon+"))*"
	elif parameter == 'Ns':
		if shuffle_type == 'off':
			Nitrogen_selection = '%.4f'%(Nitrogen_selection_unshuffled[codon])
		else:
			Nitrogen_selection = '%.4f'%(Nitrogen_selection_shuffled[codon])
		return "(math.exp("+Ns+"*"+Nitrogen_selection+"))*"
	elif parameter == 'Te':
		if shuffle_type == 'off':
			val = '%.4f'%(tAI_value[codon] - largest_value[codon_trans_standard[codon]])
		else:
			val = '%.4f'%(tAI_shuffled[codon] - largest_value[codon_trans_standard[codon]])
		return "(math.exp("+Te+"*"+val+"))*"
	else:
		print "What has gone wrong here? Not a recognised parameter %s\n" %parameter

def get_p_value(shuffle_type, real_likelihood, accuracy, Mb_in, Ns_in, Te_in): #
	counter = 1
	better_shuffled = 0
	while counter <= accuracy:
		shuf_equat = write_equation(shuffle_type)
		def equation_shuffle(x):
			Mb, Ns, Te = x
			try:
				return eval(shuf_equat)
			except OverflowError:
				return 10000000000000000000
			except ZeroDivisionError:
				return 10000000000000000000
			except ValueError:
				return 10000000000000000000
		x0 = np.asarray((Mb_in, Ns_in, Te_in))
		res_shuf = optimize.fmin(equation_shuffle, x0)
		likelihood = equation_shuffle(res_shuf)
		if likelihood <= real_likelihood:
			better_shuffled = better_shuffled + 1
		counter = counter + 1
	p_value = '%.2f'%(float(better_shuffled)/accuracy) #change this to %.3f if you change accuracy to 1000
	return p_value

def write_equation(model_type):
	global GC_unshuffled
	global Nitrogen_selection_unshuffled
	GC_unshuffled = {}
	Nitrogen_selection_unshuffled = {}
	for codon in alphabetic_codons:
		GC_unshuffled[codon] = GC_content[codon[0]] + GC_content[codon[1]] + GC_content[codon[2]]
		Nitrogen_selection_unshuffled[codon] = N_content[codon[0]] + N_content[codon[1]] + N_content[codon[2]]
	if fixed_Mb == 'moveable':
		fixed_status = 'not'
	else:
		fixed_status = fixed_Mb
	if 'shuffle' in model_type:
		global tAI_shuffled
		global GC_shuffled
		global Nitrogen_selection_shuffled
		if len(sys.argv) > 2 and "tRNAscan" in sys.argv[2]:
			tAI_shuffled = using_shuffle(tAI_value)
		GC_shuffled = using_shuffle(GC_unshuffled)
		Nitrogen_selection_shuffled = using_shuffle(Nitrogen_selection_unshuffled)
	str = "-("
	if  "Te" in model_type:
		global largest_value
		largest_value = {}
		for pro in proteins:
			largest_value[pro] = 0
		for pro in largest_value:
			for syn in proteins[pro]:
				if "shuffle_Te" in model_type:
					if tAI_shuffled[syn] > largest_value[pro]:
						largest_value[pro] = tAI_shuffled[syn]
				else:
					if tAI_value[syn] > largest_value[pro]:
						largest_value[pro] = tAI_value[syn]
	upper = {}
	for codon in alphabetic_codons:
		upper[codon] = ""
	parameters = model_type.split('_')
	number_parameters = len(parameters)
	shuffle = 'off'
	for parameter in parameters:
		number_parameters -= 1
		if parameter != 'fixed':
			if parameter == 'shuffle':
				shuffle = 'on'
			else:
				for codon in alphabetic_codons:
					if shuffle == 'off':
						additional = get_upper_equation(codon, parameter, 'Mb', 'Ns', 'Te', 'for_real', fixed_status, 'off')
					else:
						additional = get_upper_equation(codon, parameter, 'Mb', 'Ns', 'Te', 'for_shuffle', fixed_status, 'on')
					upper[codon] = upper[codon]+additional
					if number_parameters == 0:
						upper[codon] = upper[codon][:len(upper[codon])-1]
				shuffle = 'off'
	for codon in alphabetic_codons:
		prot = codon_trans_standard[codon]
		base = "("
		for syn in proteins[prot]:
			base = base+upper[syn]+" + "
		base = base[:len(base)-3]+")"
		count_codon = '%i'%codon_count[codon]
		if model_type == "Te" and prot == 'B':
			str = str+" + "
		else:
			Mb = 0.1
			Te = 0.1
			Ns = 0.1
			to_be_logged = eval(upper[codon]+"/"+base)
			if to_be_logged <=0.00000001:
				str = str+count_codon+"*(math.log(0.00000001)) + "
			elif to_be_logged >= 0.9999999:
				str = str+" + "
			else:
				str = str+count_codon+"*(math.log("+upper[codon]+"/"+base+")) + "
	str = str[:len(str)-3]+")"
	return str
	
def get_tAI_values():
	anticodons = {}
	perfect_match = {}
	for prot in proteins:
		anticodons[prot] = []
		for syn in proteins[prot]:
			perfect_match[syn] = ''
			anti = complementary_nuc[syn[2]]+complementary_nuc[syn[1]]+complementary_nuc[syn[0]]
			perfect_match[syn] = anti
			anticodons[prot].append(anti);
	anticodon_count = {}
	with open(sys.argv[2]) as tRNAscan:
		species_name = sys.argv[2]
		print species_name
		for line in tRNAscan:
			line = str(line.rstrip().strip())
			bits = line.split()
			an = bits[5]
			if "tRNA" in line or "---" in line:
				line = line
			elif "SeC" in line or an == '???':
				print "I am not including SeC finds in the analysis. Change this if you want"
			elif complementary_nuc[an[2]]+complementary_nuc[an[1]]+complementary_nuc[an[0]] in codon_trans_standard and an in anticodon_count and bits[8]>30:
				anticodon_count[an] = anticodon_count[an] + 1
			elif complementary_nuc[an[2]]+complementary_nuc[an[1]]+complementary_nuc[an[0]] in codon_trans_standard and bits[8]>30:
				anticodon_count[an] = 1
	sij_final = {}
	global tAI_value
	tAI_value = {}
	for codon in codon_trans_standard:
		tai = 0
		translated = codon_trans_standard[codon]
		for anticodon in anticodons[translated]:
			ideal = perfect_match[codon]
			if anticodon == ideal:
				sij_value = 0
			else:
				if ideal[1] == anticodon[1] and ideal[2] == anticodon[2]: #only wobble at position 3 allowed
					pair = anticodon[0]+complementary_nuc[ideal[0]]
					sij_value = sij_values[pair]
				else:
					sij_value = 0.9999
			if anticodon in anticodon_count:
				tai = tai + (1- sij_value)*anticodon_count[anticodon]
		tAI_value[codon] = tai

def using_shuffle(x):
		keys = x.keys()
		values = x.values()
		random.shuffle(values)
		return dict(zip(keys, values))

global fixed_Mb
fixed_Mb = "moveable"
if ('.fasta' in sys.argv[1]):
	model_to_use, fixed_Mb = import_sequence(sys.argv[1])
	print model_to_use
	print fixed_Mb
	if "-ind" in sys.argv:
		print "Looking at individual genes"
		if "-fix_mb" in sys.argv:
			print "using fixed Mb value from model\n"
		else:
			fixed_Mb = "moveable"
		per_gene_analysis(sys.argv[1], model_to_use)
else:
	print "The first argument must be the fasta file containing your sequences of interest\nie. python Codon_bias_analysis.py Genus_species.fasta Genus_species_tRNAscan.txt\n";
