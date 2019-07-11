import matplotlib.pyplot as plt
import glob, os
import numpy as np

#Se cuenta la cantidad de archivos
no_NAD_mirnas = len(glob.glob("Data/Adeno/NAD/miRNA/*.txt")) - 1
no_NAD_rnas = len(glob.glob("Data/Adeno/NAD/RNAseq/*.counts")) - 1

no_NSC_mirnas = len(glob.glob("Data/Squamous/NSC/miRNA/*.txt")) - 1
no_NSC_rnas = len(glob.glob("Data/Squamous/NSC/RNAseq/*.counts")) - 1

#Tejido canceroso
no_TAD_mirnas = len(glob.glob("Data/Adeno/TAD/miRNA/*.txt")) - 1
no_TAD_rnas = len(glob.glob("Data/Adeno/TAD/RNAseq/*.counts")) - 1

no_TSC_mirnas = len(glob.glob("Data/Squamous/TSC/miRNA/*.txt")) - 1
no_TSC_rnas = len(glob.glob("Data/Squamous/TSC/RNAseq/*.counts")) - 1

#Otra forma de hacerlo
#no_f = [len(files) for path, dirs, files in os.walk("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Data")]

#Desplegando informacion en pantalla
print("\nNumber of files per Tissue\n")
print("Tissue\tRNA_files\tmiRNA_files")
print("NAD\t{}\t{}".format(no_NAD_rnas,no_NAD_mirnas))
print("TAD\t{}\t{}".format(no_TAD_rnas,no_TAD_mirnas))
print("NSC\t{}\t{}".format(no_NSC_rnas,no_NSC_mirnas))
print("TSC\t{}\t{}\n".format(no_TSC_rnas,no_TSC_mirnas))
#######################
#PLOT
# Data to plot Healthy & Cancer samples
n_groups = 4
RNAseq = (no_TAD_rnas, no_TSC_rnas, no_NAD_rnas, no_NSC_rnas)
miRNASeq = (no_TAD_mirnas, no_TSC_mirnas, no_NAD_mirnas, no_NSC_mirnas)

# create plot
fig, ax = plt.subplots()
index = np.arange(n_groups)
bar_width = 0.25
opacity = 0.8

rects1 = plt.bar(index, RNAseq, bar_width,
alpha=opacity,
color='b',
label='RNAseq')

rects2 = plt.bar(index + bar_width, miRNASeq, bar_width,
alpha=opacity,
color='g',
label='miRNASeq')

plt.xlabel('Tissue Type')
plt.ylabel('no. of files')
plt.title('Number of Files of Cancer samples')
plt.xticks(index + bar_width, ('TAD', 'TSC', 'NAD', 'NSC'))
plt.legend()

plt.tight_layout()

os.chdir("Plots")
#plt.show()
plt.savefig("no_files.png")
os.chdir("..")