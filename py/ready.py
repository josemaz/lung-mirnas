import os, sys
import argparse

parser = argparse.ArgumentParser(description='Check if the previos step was done succesfully')
parser.add_argument('-s', '--step',
        help='Sample type. Ej: 1',
        required='True',
        choices=['1', '2', '3', '4'],
        default='1')
results = parser.parse_args(sys.argv[1:])

#Begining of the programm
step = results.step

#Checa que esten los archivos "index.txt" que tienen el nombre de todos los archivos descargados
if step == "1" and os.path.exists("Data/Adeno/NAD/miRNA/index.txt") and os.path.exists("Data/Adeno/NAD/RNAseq/index.txt") and os.path.exists("Data/Adeno/TAD/miRNA/index.txt") and os.path.exists("Data/Adeno/TAD/RNAseq/index.txt") and os.path.exists("Data/Squamous/NSC/RNAseq/index.txt") and os.path.exists("Data/Squamous/NSC/miRNA/index.txt") and os.path.exists("Data/Squamous/TSC/RNAseq/index.txt") and os.path.exists("Data/Squamous/TSC/miRNA/index.txt"):
	print("\nSTEP 1 Done succesfully!\n¡Everything ready for STEP 2: RNA Quality Control and Normalization!\n")
elif step == "1":
	print("\nSomething went wrong on STEP 1.\nPlease verify you have acomplished the previous step.\n")

#Checa que esten generadas la ultimas graficas despues de realizar el control de calidad
if step == "2" and os.path.exists("Plots/Adeno/RNA/QC_POST/07-PCAScore_raw.pdf") and os.path.exists("Plots/Squamous/RNA/QC_POST/07-PCAScore_raw.pdf"):
	print("\nSTEP 2 Done succesfully!\n¡Everything ready for STEP 3: miRNA Quality Control!\n")
elif step == "2":
	print("\nSomething went wrong on STEP 2.\nPlease verify you have acomplished the previous step.\n")

#Checa que esten los archivos "targets.tsv" en Adeno y en Squamous
if step == "3" and os.path.exists("Data/Adeno/rdata/miRNA/Targets.tsv") and os.path.exists("Data/Squamous/rdata/miRNA/Targets.tsv"):
	print("\nSTEP 3 Done succesfully!\n")
elif step == "3":
	print("\nSomething went wrong on STEP 3.\nPlease verify you have acomplished the previous step.\n")