import json, requests, io
import pandas as pd
import sys, argparse
import os

#Managing input
parser = argparse.ArgumentParser(description='Script to download data of lung cancer from TCGA')
parser.add_argument('-t', '--type',
        help='Sample type. Ej: NAD',
        required='True',
        choices=['NAD', 'TAD', 'NSC', 'TSC'],
        default='NAD')
results = parser.parse_args(sys.argv[1:])

#Begining of the programm
filename = results.type
fileids = pd.read_csv("Database/" + filename + ".txt", sep='\t')

cases = []
fids = []

for index, row in fileids.iterrows():

	fid = row["id"]

	with open("json/qbyfileid.json", 'r') as f:
		filters = json.load(f)
	filters['content'][0]['content']['value'] = fid

	cases_endpt = "https://api.gdc.cancer.gov/files"

	params = {
	    "filters": json.dumps(filters),
	    "fields": "cases.submitter_id",
	    "format": "TSV",
	    "size": "10000" #HACK: modificar si los casos superan los hints
	    }

	response = requests.get(cases_endpt, params = params)
	df = pd.read_csv(io.StringIO(response.text), sep='\t', header=0)
	cases.append(df.iloc[0]['cases.0.submitter_id'])
	fids.append(fid)
	print(df.iloc[0]['cases.0.submitter_id'])

df = pd.DataFrame({"case": cases, "fid": fids})

if not os.path.exists("Data"):
	os.mkdir("Data")
if not os.path.exists("Data/Adeno"):
	os.mkdir("Data/Adeno")
if not os.path.exists("Data/Squamous"):
	os.mkdir("Data/Squamous")

if filename == "NAD":
	dirname = "Adeno"
	os.mkdir("Data/Adeno/NAD")
	dirname2 = "NAD"
elif filename == "TAD":
	dirname = "Adeno"
	os.mkdir("Data/Adeno/TAD")
	dirname2 = "TAD"
elif filename == "NSC":
	dirname = "Squamous"
	os.mkdir("Data/Squamous/NSC")
	dirname2 = "NSC"
elif filename == "NSC":
	dirname = "Squamous"
	os.mkdir("Data/Squamous/TSC")
	dirname2 = "TSC"

df.to_csv("Data/" + dirname + "/" + dirname2 + "/" + filename + "-cases.tsv", sep="\t", index = False)

