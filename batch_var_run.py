'''
This will run Javi's script repeatedly for each of the samples. The
arg_list_file just needs to be a text file with the paths of each
arguments file and the script path is just the path to Javi's
script.
'''


import sys,os

arg_list_file = sys.argv[1]
script_path = sys.argv[2]

arg_list = open(arg_list_file).read().split("\n")
args = [arg for arg in arg_list if len(arg.strip()) > 0]

for arg_path in args:
	os.system("python " + script_path + " " + arg_path)
