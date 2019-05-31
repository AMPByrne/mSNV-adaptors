'''
This script will put the eight results files from Javi's script for each
sample (one per segment) together. Positions with coverage below the
user-set minimum will then be removed.It will then filter for positions 
where at least one nucleotide (not matching the reference) is present at 
a frequency above another user-defined minimum.

The parametersFile should be the arguments for this script. An example is
example.filt_args.

The aa_ref file header also needs to include the reading frame and offset
of the sequence in nt (with respect to the nucleotide reference sequence
for that gene).
'''

import sys
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

parametersFile = sys.argv[1]
for line in open(parametersFile):
	exec(line.strip())

class ref_gene:
	
	def __init__(self, name):
		self.name = name
		self.nt_seq = ""
		self.aa_seq = ""
		self.aa_off = 0
		self.aa_rf = 0
		
	def add_nt_seq(self, seq):
		self.nt_seq = seq
	
	def add_aa_seq(self, seq, ref_frame, offset):
		self.aa_seq = seq
		self.aa_rf = int(ref_frame)
		self.aa_off = int(offset)
		

def get_nt_refs(nt_ref_dir, gene_list):
	
	if not nt_ref_dir.endswith("/"):
		nt_ref_dir += "/"
	gene_paths = [nt_ref_dir + gene for gene in gene_list]
	gene_seqs = ["".join(open(path).read().split("\n")[1:-1]) for path in gene_paths]
	ref_dict = {gene: seq for gene, seq in zip(gene_list, gene_seqs)}
	
	return ref_dict
	
def get_aa_refs(aa_ref_dir, gene_list):
	
	if not aa_ref_dir.endswith("/"):
		aa_ref_dir += "/"
	aa_paths = [aa_ref_dir + gene for gene in gene_list]
	aa_data = [open(path).read().split("\n")[:-1] for path in aa_paths]
	aa_info = [data[0].split()[-1].split("_")[1:] for data in aa_data]
	aa_seqs = ["".join(open(path).read().split("\n")[1:-1]) for path in aa_paths]
	comb = [[seq] + info for seq, info in zip(aa_seqs, aa_info)]
	ref_dict = {gene: comb_info for gene, comb_info in zip(gene_list, comb)}
	
	return ref_dict
	
def init_ref_inst(name, nt_ref_dict, aa_ref_dict):
	
	aa_seq = aa_ref_dict[name][0]
	aa_rf = aa_ref_dict[name][1]
	aa_off = aa_ref_dict[name][2]
	
	class_inst = ref_gene(name)
	class_inst.add_nt_seq(nt_ref_dict[name])
	class_inst.add_aa_seq(aa_seq, aa_rf, aa_off)
	
	return class_inst

def get_results(res_path, gene_list):
	
	res_paths = [res_path.replace("###",gene) for gene in gene_list]
	res_list = [pd.read_csv(path) for path in res_paths]
	comb_res = pd.concat(res_list)
	
	return comb_res

def fill_df(row, ref_dict, min_freq):
	
	#Filling codon positions and ref AAs
	gene = row["gene"]
	pos = row["pos"] - 1
	rf = ref_dict[gene].aa_rf - 1
	offset = ref_dict[gene].aa_off
	codon_pos = ((pos - rf) // 3) - offset
	if codon_pos < 0:
		codon_pos = np.nan
	row["codon"] = codon_pos + 1
	try:
		ref_AA = ref_dict[gene].aa_seq[int(codon_pos)]
	except:
		ref_AA = np.nan
	row["ref_AA"] = ref_AA
	
	#Filling ref nts
	pos = row["pos"] - 1
	ref_nt_seq = ref_dict[gene].nt_seq
	ref_base = ref_nt_seq[pos]
	row["ref_nt"] = ref_base
	
	#Ordering variant nts, filling nt freqs,
	#translating variant codons and testing if variants
	#are synonymous
	cov = row["cov"]
	bases = ["A","C","G","T"]
	if cov > 0:
		base_freqs = {base: row[base]/cov for base in bases}
	else:
		base_freqs = {base: 0 for base in bases}
	sorted_base_freqs = sorted(base_freqs.items(), reverse=True, key=lambda x: x[1])
	pos_in_codon = (pos - rf) % 3
	codon = list(ref_nt_seq[pos - pos_in_codon : pos + (3 - pos_in_codon)])
	for num in range(1,5):
		i = num - 1
		base = sorted_base_freqs[i][0]
		freq = round(sorted_base_freqs[i][1], 3)
		row["nt_" + str(num)] = base
		row["freq_" + str(num)] = freq
		try:
			codon[pos_in_codon] = base
		except:
			None
		s_codon = Seq("".join(codon), IUPAC.ambiguous_dna)
		if len(s_codon) == 3:
			var_AA = str(s_codon.translate())
			if float(freq) >= min_freq and type(ref_AA) == str:
				if var_AA == ref_AA:
					row["syn_" + str(num)] = "syn"
				else:
					row["syn_" + str(num)] = "non-syn"
		else:
			var_AA = np.nan
		row["AA_" + str(num)] = var_AA
	
	return row

#Reading in reference sequences
genes = ["HA", "MA", "NA", "NP", "NS", "PA", "PB1", "PB2"]
nt_refs = get_nt_refs(nt_ref_path, genes)
aa_refs = get_aa_refs(aa_ref_path, genes)
refs = {gene: init_ref_inst(gene, nt_refs, aa_refs) for gene in genes}

#Reading in inputs
results = get_results(results_path, genes)
#this is because pandas interprets the gene name "NA" as empty
results["seq"] = results["seq"].fillna("NA")
results = results[results["afterFilt"] >= min_cov]

#Creating new dataframe for results
cols = ["gene","pos","codon","cov","ref_nt","ref_AA","",
		"freq_1","nt_1","AA_1","syn_1","",
		"freq_2","nt_2","AA_2","syn_2","",
		"freq_3","nt_3","AA_3","syn_3","",
		"freq_4","nt_4","AA_4","syn_4","",
		"A","C","G","T"]
df = pd.DataFrame(columns = cols)
df["gene"] = results["seq"]
df["pos"] = results["pos"]
df["cov"] = results["afterFilt"]
bases = ["A", "C", "G", "T"]
for base in bases:
	df[base] = results[base]
df = df.apply(fill_df, axis = 1, ref_dict = refs, min_freq = min_var_freq)
for base in bases:
	del df[base]
df = df[(df["ref_nt"] != df["nt_1"]) | (df["freq_2"] >= min_var_freq)]


#Saving outputs
sample_name = results_path.split("/")[-1].split("_")[0]
if not output_path.endswith("/"):
	output_path += "/"
all_filtered_res_path = output_path + sample_name + "_filtered_variants.csv"
df.to_csv(all_filtered_res_path)

