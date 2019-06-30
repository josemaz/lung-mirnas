#libraries
import pandas as pd
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt

#abrimos los archivos
file1 = pd.read_csv("Data/NAD-mirna.tsv", sep="\t")
file2 = pd.read_csv("Data/NSC-mirna.tsv", sep="\t")
file3 = pd.read_csv("Data/TAD-mirna.tsv", sep="\t")
file4 = pd.read_csv("Data/TSC-mirna.tsv", sep="\t")

all_mirnas = []
all_mirnas.append(list(file1["mirna_count"]))
all_mirnas.append(list(file2["mirna_count"]))
all_mirnas.append(list(file3["mirna_count"]))
all_mirnas.append(list(file4["mirna_count"]))

f_names = ["NAD", "NSC", "TAD", "TSC"]
cero_mirnas = []
one_mirnas = []
two_mirnas = []
three_mirnas = []

#We make the plot for every group
i = 0
for (element, f_name) in zip(all_mirnas, f_names):
	print(f_name)
	cuentas = {number: element.count(number) for number in element}
	print(cuentas)

	#Guardamos en listas la cantidad de mirnas de cada tipo que tiene cada archivo
	#Se va a usar en el resumen gral
	if len(cuentas) == 2:
		cero_mirnas.append(0)	
		one_mirnas.append(cuentas.values()[0])	
		two_mirnas.append(cuentas.values()[1])
	else:
		cero_mirnas.append(cuentas.values()[0])	
		one_mirnas.append(cuentas.values()[1])
		two_mirnas.append(cuentas.values()[2])
	if len(cuentas) == 4:
		three_mirnas.append(cuentas.values()[3])
	else:
		three_mirnas.append(0)	
	###########
	height = cuentas.values()
	if len(height) == 4:
		bars = ('0','1','2','3')
	elif len(height) == 3:
		bars = ('0','1','2')
	elif len(height) == 2:
		bars = ('1','2')

	y_pos = np.arange(len(bars))
	colors = ["black","blue","green","red"]

	plt.bar(bars, height,color = colors[i], align='center')
	plt.title('Frecuency of miRNAS in ' + f_names[i])
	plt.xticks(y_pos, bars)
	plt.xlabel('no. miRNAs')
	plt.ylabel('Frequency')
	plt.show()
	i+=1

#This plot is a Resume for all the previous information
n_groups = 4

fig, ax = plt.subplots()
index = np.arange(n_groups)
bar_width = 0.20
opacity = 0.8

rects1 = plt.bar(index, cero_mirnas, bar_width,
alpha=opacity,
color='k',
label='0 miRNAs')

rects2 = plt.bar(index + bar_width, one_mirnas, bar_width,
alpha=opacity,
color='b',
label='1 miRNA')

rects3 = plt.bar(index + bar_width + bar_width, two_mirnas, bar_width,
alpha=opacity,
color='g',
label='2 miRNAs')

rects4 = plt.bar(index + bar_width + bar_width + bar_width, three_mirnas, bar_width,
alpha=opacity,
color='r',
label='3 miRNAs')

plt.xlabel('Tissue Type')
plt.ylabel('no. of miRNAs')
plt.title('Number of miRNAs of Healthy & Cancerous samples')
plt.xticks(index + bar_width, ('NAD', 'NSC', 'TAD', 'TSC'))
plt.legend()

plt.tight_layout()
plt.show()