#libraries
import pandas as pd
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import os

#abrimos los archivos
file1 = pd.read_csv("Data/NAD-mirna.tsv", sep="\t")
file2 = pd.read_csv("Data/NSC-mirna.tsv", sep="\t")
file3 = pd.read_csv("Data/TAD-mirna.tsv", sep="\t")
file4 = pd.read_csv("Data/TSC-mirna.tsv", sep="\t")

all_mirnas = []
all_mirnas.append(list(file1["mirna_count"])) #lista con las cuentas
all_mirnas.append(list(file2["mirna_count"]))
all_mirnas.append(list(file3["mirna_count"]))
all_mirnas.append(list(file4["mirna_count"]))

f_names = ["NAD", "NSC", "TAD", "TSC"]
cero_mirnas = []
one_mirnas = []
two_mirnas = []
three_mirnas = []
number_mirnas = [cero_mirnas,one_mirnas,two_mirnas,three_mirnas]

os.mkdir("Graphs")
os.chdir("Graphs")

#We make the plot for every group
i = 0
for (element, f_name) in zip(all_mirnas, f_names):
	cuentas = {number: element.count(number) for number in element}
	
	listita = cuentas.items() #Tiene una lista de tuplas [(0,veces),(1,veces)]

	for tupla in listita:
		if tupla[0] == 0:
			cero_mirnas.append(tupla[1])
		if tupla[0] == 1:
			one_mirnas.append(tupla[1])
		if tupla[0] == 2:
			two_mirnas.append(tupla[1])
		if tupla[0] == 3:
			three_mirnas.append(tupla[1])
	
	#No en todos los casos va haber con 0,1,2 o 3 miRNAs, entonces buscas esos casos y pones un cer en esa posicion
	opciones = [0,1,2,3]  
	for i in range(len(opciones)):
		if opciones[i] not in list(cuentas.keys()):
			number_mirnas[i].insert(opciones[i],0)

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
	#plt.show()
	#plt.savefig("no_mirnas_" + f_names[i] + ".png")
	i+=1

#Desplegamos en la pantalla los resultados
print("\nmiRNAs p/RNA   NAD   TAD   NSC   TSC")
print("     0          {}    {}    {}    {}".format(cero_mirnas[0],cero_mirnas[1],cero_mirnas[2],cero_mirnas[3]))
print("     1          {}    {}    {}    {}".format(one_mirnas[0],one_mirnas[1],one_mirnas[2],one_mirnas[3]))
print("     2          {}    {}    {}    {}".format(two_mirnas[0],two_mirnas[1],two_mirnas[2],two_mirnas[3]))
print("     3          {}    {}    {}    {}\n".format(three_mirnas[0],three_mirnas[1],three_mirnas[2],three_mirnas[3]))


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
#plt.show()
plt.savefig("no_mirnas_ALL.png")

os.chdir("..")