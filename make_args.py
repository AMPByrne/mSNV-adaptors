'''
This helps make arguments files for each run of Javi's script - 
needs to be done for batch_var_run.py to be used.

The sample_file needs to be a text file containing the path to each of
the reads files. Only put one line per sample, and replace R1/R2 with ###
and the script will put the R1 and R2 paths in the right place in the arguments
file.

The ref_file needs to be another text file containing the paths of each
of the reference sequences the reads will be mapped to. In this case
it's basically going to be the paths of eight separate fasta files, each
containing the reference sequence for one of the segments.

The args_dir should be an empty directory where the script should put all
the new arguments files.

The output_dir is the directory where you want Javi's script to put its
results (probably on mounted storage as some of the files get quite large).

Finally, this script will then create a text file containing the path
to each of the new arguments files for use by batch_var_run.py. The 
arg_list_fname is what you want this file to be called.
'''

import sys

sample_file = sys.argv[1]
ref_file = sys.argv[2]
args_dir = sys.argv[3]
output_dir = sys.argv[4]
arg_list_fname = sys.argv[5]

samples = open(sample_file).read().split("\n")[:-1]
refs = open(ref_file).read().split("\n")[:-1]
arg_list = []

pairs = [(samp,ref) for samp in samples for ref in refs]
for (samp,ref) in pairs:
	samp_name = samp.split("/")[-1].split("_")[0]
	ref_name = ref.split("/")[-1]
	if not args_dir.endswith("/"):
		args_dir += "/"
	arg_fname = args_dir + samp_name + "_" + ref_name + ".args"
	arg_list.append(arg_fname)
	R1 = samp.replace("###","R1")
	R2 = samp.replace("###","R2")
	with open(arg_fname, "w") as fh:
		fh.write("ref_sequence=\"" + ref + "\"\n")
		fh.write("fastq_R1=\"" + R1 + "\"\n")
		fh.write("fastq_R2=\"" + R2 + "\"\n")
		fh.write("call_quality_threshold=30\n")
		fh.write("map_quality_threshold=40\n")
		fh.write("software_path=\"/home/p946572/python/Javi_stuff\"\n")
		fh.write("results_path=\"" + output_dir + "\"\n")
		fh.write("trimmomatric=\"trimmomatic\"\n")
		fh.write("vcfutils=\"vcfutils.pl\"\n")

with open(arg_list_fname,"w") as fh:
	fh.write("\n".join(arg_list))
