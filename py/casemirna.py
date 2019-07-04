import json, requests, io
import pandas as pd
import numpy as np
import sys, argparse

#Managing input
parser = argparse.ArgumentParser(description='Script to download data of lung cancer from TCGA')
parser.add_argument('-t', '--type',
        help='Sample type. Ej: NAD',
        required='True',
        choices=['NAD', 'TAD', 'NSC', 'TSC'],
        default='NAD')
results = parser.parse_args(sys.argv[1:])

filename = results.type

if filename == "NAD":
	dirname = "Adeno"
	dirname2 = "NAD"
elif filename == "TAD":
	dirname = "Adeno"
	dirname2 = "TAD"
elif filename == "NSC":
	dirname = "Squamous"
	dirname2 = "NSC"
elif filename == "NSC":
	dirname = "Squamous"
	dirname2 = "TSC"

#Begining of the programm
cases = pd.read_csv("Data/" + dirname + "/" + dirname2 + "/" + filename + "-cases.tsv", sep='\t')
mirna_fid = []
mirna_fname = []
mirna_count = []
for index, row in cases.iterrows():

	print(row['case'])

	with open("json/qbyMIRNA.json", 'r') as f:
		filters = json.load(f)
	filters['content'][0]['content']['value'] = row['case']

	cases_endpt = "https://api.gdc.cancer.gov/files"

	params = {
	    "filters": json.dumps(filters),
	    "fields": "file_name,data_format,file_id",
	    "format": "TSV",
	    "size": "10000" #HACK: modificar si los casos superan los hints
	    }
	response = requests.get(cases_endpt, params = params)
	try:
		df = pd.read_csv(io.StringIO(response.text), sep='\t', header=0)
		mirna_count.append(df.shape[0])
		mirna_fid.append(df.loc[0, "file_id"])
		mirna_fname.append(df.loc[0, 'file_name'])
	except:
		df = np.nan
		mirna_count.append(0)
		mirna_fid.append(df)
		mirna_fname.append(df)

cases['mirna_count'] = mirna_count
cases['mirna_fname'] = mirna_fname
cases['mirna_fid'] = mirna_fid
cases.rename(columns={'fid':'rnaseq_fid'}, inplace = True)

cases.to_csv("Data/" + dirname + "/" + dirname2 + "/" + filename + "-mirna.tsv", sep="\t", index = False)

