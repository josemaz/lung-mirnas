import json, requests, io
import pandas as pd


filename = sys.argv[1]
fileids = pd.read_csv("../Database/" + filename + ".txt", sep='\t')

cases = []
fids = []

for index, row in fileids.iterrows():

	fid = row["id"]

	with open("../json/qbyfileid.json", 'r') as f:
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
df.to_csv("../Data/" + filename + "-cases.tsv", sep="\t", index = False)

