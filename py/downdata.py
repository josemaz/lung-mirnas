import pandas as pd
import sys, os, shutil
import Util
from colored import fore, back, style

dirname = sys.argv[1]

db = pd.read_csv("../Database/" + dirname + "-mirna.tsv", sep='\t')
db = db[db.mirna_count != 0]

if os.path.exists(dirname):
	shutil.rmtree(dirname)
os.mkdir(dirname)
os.chdir(dirname)

if os.path.exists("RNAseq"):
	shutil.rmtree('RNAseq')
os.mkdir("RNAseq")
if os.path.exists("miRNA"):
	shutil.rmtree('miRNA')
os.mkdir("miRNA")

print(fore.LIGHT_BLUE + style.BOLD + "Downloading RNAseq data..." + style.RESET)
os.chdir("RNAseq")
for i,row in db.iterrows():
	Util.download(row['rnaseq_fid'])
os.chdir("..")

print(fore.LIGHT_BLUE + style.BOLD + "Downloading miRNA data..." + style.RESET)
os.chdir("miRNA")
for i,row in db.iterrows():
	Util.download(row['mirna_fid'])
os.chdir("..")
