# lung-mirnas
Proyecto para procesar RNAseq y miRNAs para cáncer de pulmón

CONTENIDO DE CARPETAS

**Database (12 files)**   
   Contiene de NAD,TAD,NSC & TSC, 3 archivos de c/u. El ".*.txt" es un archivo manifest, el ".*-cases.tsv" contiene el caseID y el fid, el ".*-mirna.tsv" contiene caseID, rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid

**R (2 files)**   
   Tiene script para el manejo de RNASeq raw counts

**json (2 files)**   
   Tiene un archivo con el que se hizo el query a TCGA y obtener los CaseID a partir de los manifest;  Tiene otro archivo para hacer un query a TCGA y obtener rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid a partir de los CaseID

**py (6 files)**   
*script 1: Util.py   
   Libreria creada para el download de informacion   
*script 2: casemirna.py   
   Script para obtener rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid a partir de los CaseID
*script 3: downdata.py   
   Script para bajar archivos de cuentas de RNASeq y de miRNASeq y dentro de crapetas, ademas crear un archivo llamado "index.txt" con un listado de tosos los archivos descargados.   
*script 4: file_number.txt   
   Script para saber el numero de archivos de raw counts se tienen tanto de RNASeq como de miRNASeq y hacer graficas
*script 5: getcases.py   
   Script para obtener los caseID a partir de archivos manifest   
*script 6: getmirnasCounts.py   
   Script para graficar el numero de miRNAs de NAD,TAD,NSC y de TSC



MANUAL PARA DESCARGAR DATOS

Se partio de archivos manifest.   
De estos se ejecuta el "script 5: getcases.py" de python para que a partir de de esos archivos se obtenga como output archivos ".*-cases.tsv" con el caseID y el fid.

Se ejecuta el "script 5: getcases.py" para que a partir de los archivos ".*-cases.tsv" se obtengan archivos a ".*-mirna.tsv" "" donde haya caseID, rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid.
Se ejecuta el "script 2: casemirna.py" para que a partir de archivos ".*-mirna.tsv" se obtengan archivos ".*-mirna.tsv" con caseID, rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid. 

Una vez con los archivos anteriores se ejecuta el "script 4: file_number.txt" & el "script 6: getmirnasCounts.py" para obtener graficas que ilustren los mirnas y la cantidad dearchivos.

Se ejecuta el "script 3: downdata.py" para descargar de TCGA los archivos con las raw counts de RNASeq y de miRNASeq y generar el archivo "index.txt" con un listado de los archivos descargados.

Se ejecuta el script de R "Rectify_rnaSeq.R" para manejar datos.
