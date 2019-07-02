import requests, json, re
import gzip, shutil, os
import pandas as pd
import io, sys
import Util, argparse
#from colored import fore, back, style

#Usefull function
def do_file (dir_name):
	if os.path.exists(dir_name):
		shutil.rmtree(dir_name)
	os.mkdir(dir_name)

#Managing input
parser = argparse.ArgumentParser(description='Script to download data of lung cancer from TCGA')
parser.add_argument('-t', '--type',
        help='Sample type. Ej: NAD',
        required='True',
        choices=['NAD', 'TAD', 'NSC', 'TSC'],
        default='NAD')
results = parser.parse_args(sys.argv[1:])

#Begining of the programm
dirname = results.type
db = pd.read_csv("Data/" + dirname + "-mirna.tsv", sep='\t')
db = db[db.mirna_count != 0]

os.chdir('Data/')
do_file(dirname)
os.chdir(dirname)

do_file('RNAseq')

do_file('miRNA')

#print(fore.LIGHT_BLUE + style.BOLD + "Downloading RNAseq data..." + style.RESET)
print("\nDownloading RNAseq data...")
os.chdir("RNAseq")

with open("index.txt", mode='w+') as rna_file: 
	for i,row in db.iterrows():		
		print(i)
		file_name = Util.download(row['rnaseq_fid'])
		case_ID = row["case"]
		rna_file.write(case_ID + "\t" + file_name + '\n')
os.chdir("..")

#print(fore.LIGHT_BLUE + style.BOLD + "Downloading miRNA data..." + style.RESET)
print("\nDownloading miRNA data...")
os.chdir("miRNA")

with open("index.txt", mode='w+') as mirna_file:
	for i,row in db.iterrows():
		print(i)
		file_name = Util.download(row['mirna_fid'])
		case_ID = row["case"]
		mirna_file.write(case_ID + "\t" + file_name + '\n')
os.chdir("..")
