# lung-mirnas
Proyecto para procesar RNAseq y miRNAs para cáncer de pulmón

## Contenido deCarpetas

- **Database (4 files)**    
   Esta carpeta es **fija** y contine los archivos originales de cancer de pulmon 
   (NAD.txt, TAD.txt, NSC.txt, TCS,TXT).

- **Data (8 files, 4 directories)**  
   Contiene de NAD,TAD,NSC & TSC,  archivos de c/u. El ".*-cases.tsv" contiene el caseID y el fid, el ".*-mirna.tsv" contiene caseID, rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid

- **R (2 files)**  
   Tiene script para el manejo de RNASeq raw counts

- **json (2 files)**  
   Tiene un archivo con el que se hizo el query a TCGA y obtener los CaseID a partir de los manifest;  Tiene otro archivo para hacer un query a TCGA y obtener rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid a partir de los CaseID

- **py (6 files)**  
   - Util.py: Libreria creada para el download de informacion
   - casemirna.py: Script para obtener rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid a partir de los CaseID
   - downdata.py: Script para bajar archivos de cuentas de RNASeq y de miRNASeq y dentro de crapetas, ademas crear un archivo llamado "index.txt" con un listado de tosos los archivos descargados.   
   - file_number.txt: Script para saber el numero de archivos de raw counts se tienen tanto de RNASeq como de miRNASeq y hacer graficas
   -getcases.py: Script para obtener los caseID a partir de archivos manifest
   -getmirnasCounts.py: Script para graficar el numero de miRNAs de NAD,TAD,NSC y de TSC

# Instrucciones para procesar los datos

## Prerequistos
 - Python (3.7.3)
 - R (3.6.0)
 - librerias (matplotlib.pyplot, numpy, glob, pandas, json, requests, re, gzip, shutil, )

## Descarga de los datos (python)
   1. `python getcases.py NAD` [Obtencion de Case ID para poder bajar datos]
   2. `python getcases.py TAD`
   3. `python getcases.py NSC`
   4. `python getcases.py TSC`
   ...
   5. `python casemirna.py NAD` [Obtencion de mirnafile_ID para bajar datos]
   6. `python casemirna.py TAD`
   7. `python casemirna.py NSC`
   8. `python casemirna.py TSC`
   ...
   9. `python getmirnasCounts.py` [Grafica de la cantidad de miRNAs]
   ...
   10. `python downdata.py TSC` [Breve descripcion]
   11. `python downdata.py NAD`
   12. `python downdata.py TAD`
   13. `python downdata.py NSC`
   ...
   13. `python file_number.py` [Grafica de la cantidad de archivos]
   

## Control de calidad y normalizacion (R)
   1. `Rscript 01-Rectify_rnaSeq.R` [Breve descripcion]
   2. `Rscript 02-PRE-QC.R` [Breve descripcion]

## Descarga
 
 - __getcases.py__
 
 
Se partio de archivos manifest.   
De estos se ejecuta el "script 5: getcases.py" de python para que a partir de de esos archivos se obtenga como output archivos ".*-cases.tsv" con el caseID y el fid.

Se ejecuta el "script 5: getcases.py" para que a partir de los archivos ".*-cases.tsv" se obtengan archivos a ".*-mirna.tsv" "" donde haya caseID, rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid.
Se ejecuta el "script 2: casemirna.py" para que a partir de archivos ".*-mirna.tsv" se obtengan archivos ".*-mirna.tsv" con caseID, rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid. 

Una vez con los archivos anteriores se ejecuta el "script 4: file_number.txt" & el "script 6: getmirnasCounts.py" para obtener graficas que ilustren los mirnas y la cantidad dearchivos.

Se ejecuta el "script 3: downdata.py" para descargar de TCGA los archivos con las raw counts de RNASeq y de miRNASeq y generar el archivo "index.txt" con un listado de los archivos descargados.

Se ejecuta el script de R "Rectify_rnaSeq.R" para manejar datos.

