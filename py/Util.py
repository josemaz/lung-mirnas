import requests, json, re
import gzip, shutil, os
import pandas as pd
import io, sys

#####################################################
### Funcion para descargar un listado
def download(fid):
#file_id = "b658d635-258a-4f6f-8377-767a43771fe4"
#file_id = "26562f3a-6491-477c-888b-8cd9e994233a"

	data_endpt = "https://api.gdc.cancer.gov/data/{}".format(fid)
	response = requests.get(data_endpt, headers = {"Content-Type": "application/json"})
	response_head_cd = response.headers["Content-Disposition"]
	file_name = re.findall("filename=(.+)", response_head_cd)[0]
	with open(file_name, "wb") as output_file:
    		output_file.write(response.content)
	fext = os.path.splitext(file_name)[1]
	print(file_name)
	if fext == ".gz":
		cmd = "gunzip " + file_name
		os.system(cmd)

	return(file_name)