# lung-mirnas
Proyecto para procesar RNAseq y miRNAs para cáncer de pulmón

## Contenido deCarpetas

- **Database (4 files)**    
   Esta carpeta es **fija** y contine los archivos originales de cancer de pulmón (NAD.txt, TAD.txt, NSC.txt, TCS,TXT).

- **Data (8 files, 4 directories)**  
   sta carpeta es **generada** y contiene de NAD,TAD,NSC & TSC archivos de c/u. El ".*-cases.tsv" contiene el caseID y el fid, el ".*-mirna.tsv" contiene caseID, rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid  
   
- **Graphs (2 files)**  

- **json (2 files)**  
   Esta carpeta es **fija** y tiene un archivo con el que se hizo el query a TCGA y obtener los CaseID a partir de los manifest;  Tiene otro archivo para hacer un query a TCGA y obtener rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid a partir de los CaseID
   - qbyfileid.json: Script para 
   - qbyMIRNA.json: Script para 
   
- **pipeline (1 file)**  

- **py (6 files)**  
   - Util.py: Libreria creada para el download de informacion
   - casemirna.py: Script para obtener rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid a partir de los CaseID
   - downdata.py: Script para bajar archivos de cuentas de RNASeq y de miRNASeq y dentro de crapetas, ademas crear un archivo llamado "index.txt" con un listado de tosos los archivos descargados.   
   - file_number.txt: Script para saber el numero de archivos de raw counts se tienen tanto de RNASeq como de miRNASeq y hacer graficas  
   - getcases.py: Script para obtener los caseID a partir de archivos manifest
   - getmirnasCounts.py: Script para graficar el numero de miRNAs de NAD,TAD,NSC y de TSC
  
- **R (4 files)**  
   - 01-Rectify_rnaSeq.R: Script para 
   - 02-PRE-QC.R: Script para 
   - 03-NORM.R: Script para 
   - 04-POST-QC.R: Script para 

   

# Instrucciones para procesar los datos

## Prerequistos
 - Python (3.7.3)
 - librerias python (matplotlib.pyplot, numpy, glob, pandas, json, requests, re, gzip, shutil)
 - R (3.6.0)
 - librerias R (BiocParallel, parallel, NOISeq, EDASeq)

## Descarga de los datos (python)
   1. `bash sh/downdata.sh` [Descarga de datos de RNASeq y miRNASeq; Gráfica de miRNAs y de Archivos]

## Control de calidad y normalizacion (R)
   1. `Rscript R/01-Rectify_rnaSeq.R` [Preparacion de archivos]
   2. `Rscript R/02-PRE-QC.R` [Control de calidad, graficas de los datos obtenidos]
   3. `Rscript R/03-NORM.R` [Normalización]


 
 

