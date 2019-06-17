import pandas as pd
import matplotlib.pyplot as plt

file1 = pd.read_csv("../Data/NAD-mirna.tsv", sep="\t")
file2 = pd.read_csv("../Data/NSC-mirna.tsv", sep="\t")
file3 = pd.read_csv("../Data/TAD-mirna.tsv", sep="\t")
file4 = pd.read_csv("../Data/TSC-mirna.tsv", sep="\t")

all_mirnas = []
all_mirnas.append(list(file1["mirna_count"]))
all_mirnas.append(list(file2["mirna_count"]))
all_mirnas.append(list(file3["mirna_count"]))
all_mirnas.append(list(file4["mirna_count"]))

f_names = ["NAD", "NSC", "TAD", "TSC"]

for (element, f_name) in zip(all_mirnas, f_names):
	print(f_name)
	cuentas = {number: element.count(number) for number in element}
	print(cuentas)
	n, bins, patches = plt.hist(x= element, bins='auto', color='#0504aa',
                            alpha=0.75, rwidth=0.85)
	plt.xlabel("mirnas counts")
	plt.ylabel("Frequency")
	plt.title("Frequency of number of mirnas.")
	plt.grid(axis="y", alpha =0.75)
	plt.show()

