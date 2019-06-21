import matplotlib.pyplot as plt
import glob
import numpy as np

#Se cuenta la cantidad de archivos
no_NAD_mirnas = len(glob.glob("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/NAD/miRNA/*.txt"))
no_NAD_rnas = len(glob.glob("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/NAD/RNAseq/*.counts"))

no_NSC_mirnas = len(glob.glob("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/NSC/miRNA/*.txt"))
no_NSC_rnas = len(glob.glob("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/NSC/RNAseq/*.counts"))

#Tejido canceroso
no_TAD_mirnas = len(glob.glob("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/TAD/miRNA/*.txt"))
no_TAD_rnas = len(glob.glob("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/TAD/RNAseq/*.counts"))

no_TSC_mirnas = len(glob.glob("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/TSC/miRNA/*.txt"))
no_TSC_rnas = len(glob.glob("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/TSC/RNAseq/*.counts"))

#Otra forma de hacerlo
#no_f = [len(files) for path, dirs, files in os.walk("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database")]

#######################
#PLOT
# Data to plot Healthy samples
n_groups = 2
RNAseq = (no_NAD_rnas, no_NSC_rnas)
miRNASeq = (no_NAD_mirnas, no_NSC_mirnas)

# create plot
fig, ax = plt.subplots()
index = np.arange(n_groups)
bar_width = 0.35
opacity = 0.8

rects1 = plt.bar(index, RNAseq, bar_width,
alpha=opacity,
color='b',
label='RNAseq')

rects2 = plt.bar(index + bar_width, miRNASeq, bar_width,
alpha=opacity,
color='g',
label='miRNASeq')

plt.xlabel(' Type')
plt.ylabel('no. of files')
plt.title('Number of Files of Healthy samples')
plt.xticks(index + bar_width, ('NAD', 'NSC'))
plt.legend()

plt.tight_layout()
plt.show()


# Data to plot for Cancer samples
n_groups = 2
RNAseq = (no_TAD_rnas, no_TSC_rnas)
miRNASeq = (no_TAD_mirnas, no_TSC_mirnas)

# create plot
fig, ax = plt.subplots()
index = np.arange(n_groups)
bar_width = 0.35
opacity = 0.8

rects1 = plt.bar(index, RNAseq, bar_width,
alpha=opacity,
color='b',
label='RNAseq')

rects2 = plt.bar(index + bar_width, miRNASeq, bar_width,
alpha=opacity,
color='g',
label='miRNASeq')

plt.xlabel('Cancer Type')
plt.ylabel('no. of files')
plt.title('Number of Files of Cancer samples')
plt.xticks(index + bar_width, ('TAD', 'TSC'))
plt.legend()

plt.tight_layout()
plt.show()


