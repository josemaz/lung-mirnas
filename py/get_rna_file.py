import json, requests, io
import pandas as pd
import numpy as np

#/home/kevinml/Documentos/INMEGEN/lung-mirnas
cases = pd.read_csv("../Data/NAD-mirna.tsv", sep='\t')
rna_fname =[]

for index, row in cases.iterrows():

	print(row['case'])

	with open("../json/get_rna_fname.json", 'r') as f:
		filters = json.load(f)
	filters['content'][0]['content']['value'] = row['case']

	cases_endpt = "https://api.gdc.cancer.gov/files"

	params = {
	    "filters": json.dumps(filters),
	    "fields": "file_name,data_format",
	    "format": "TSV",
	    "size": "10000" #HACK: modificar si los casos superan los hints
	    }
	response = requests.get(cases_endpt, params = params)
	try:
		df = pd.read_csv(io.StringIO(response.text), sep='\t', header=0)
		rna_fname.append(df.loc[0, 'file_name'])
	except:
		df = np.nan
		rna_fname.append(df)

cases['rna_fname'] = rna_fname

cases.to_csv("../Data/NAD-mirna_and_rna.tsv", sep="\t", index = False)