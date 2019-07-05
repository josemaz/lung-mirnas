# lung-mirnas
Proyecto para procesar RNAseq y miRNAs para cáncer de pulmón.

## Contenido de Carpetas

- **Database (4 files)**    
   Esta carpeta es **fija** y contine los archivos originales (manifest) de los diferentes tipos de tejido de cancer de pulmón (NAD.txt, TAD.txt, NSC.txt, TSC,txt).

- **Data (2 directories)**  
   Carpeta **generada**. Contiene 2 directorios: "Adeno" & "Squamous", que contienen toda la informacion porducida de estos dos  grupos de tejido. 
     
   Hay 2 directorios (tejido Normal y Tumoral) en cada uno de estos directorios; cada uno contiene un archivo de cuentas de RNASeq y otro de miRNASeq.
   
- **Graphs (2 files)**  
   Carpeta **generada** con 2 archivos ".png".
   - no_mirnas_ALL.png: Gráfica de cantidad de miRNAs por RNA de los 4 tipos de tejido.
   - no_files.png: Gráfica de cantidad de archivos descargados de RNASeq y de miRNASeq de los 4 tipos de tejido.

- **json (2 files)**  
   Carpeta **fija** con 2 scripts.
   - qbyfileid.json: Script para hacer el query a "The Cancer Genome Atlas" (TCGA) y obtener los CaseID a partir de los archivos manifest.
   - qbyMIRNA.json: Script que ayudara a obtener los nombres de los archivos a descargar a partir de CaseID.
   
- **pipeline (1 file)**  
   Carpeta **fija** con 1 archivo.
   - biomart-20181212.txt: Archivo con la anotación de todos los genes de *Homo sapiens*, a partir de donde se mapearan los genes  obtenidos de los archivos de RNASeq. 

- **py (6 files)**  
   Carpeta **fija** con 6 scripts.
   - Util.py: Libreria creada para la descarga de información.
   - casemirna.py: Script para obtener rnaseq_fid, cantidad de mirnas, mirna_fname y el mirna_fid a partir de los CaseID
   - downdata.py: Script para bajar archivos de cuentas de RNASeq y de miRNASeq y dentro de crapetas, ademas crear un archivo llamado "index.txt" con un listado de tosos los archivos descargados.   
   - file_number.txt: Script para saber el numero de archivos de raw counts se tienen tanto de RNASeq como de miRNASeq y hacer graficas  
   - getcases.py: Script para obtener los caseID a partir de archivos manifest
   - getmirnasCounts.py: Script para graficar el numero de miRNAs de NAD,TAD,NSC y de TSC
  
- **R (4 files)**  
   Carpeta **fija** con 4 scripts.
   - 01-Rectify_rnaSeq.R: Script para checar si las muestras tienen el mismo tamaño, checar si los genes mapean a posiciones, cambiar la anotación, remover aquellos genes mapeados a cromosomas no convencionales, remover aquellos genes que no tienen un símbolo y salvar la información!
   - 02-PRE-QC.R: Script que hace un pre-control de calidad a los datos y genera plots donde se muestra el biotipo tipo de genes obtenidos y valores de expresión.
   - 03-NORM.R: Script para realizar la normalización de la información.
   - 04-POST-QC.R: Script que hace un post-control de calidad a los datos.

   

# Instrucciones para procesar los datos

## Prerrequistos
 - Python (3.7.3)
 - librerias python (matplotlib.pyplot, numpy, glob, pandas, json, requests, re, gzip, shutil)
 - R (3.6.0)
 - librerias R (BiocParallel, parallel, NOISeq, EDASeq)

## Descarga de los datos (se usan los scripts de python)
   1. `bash sh/downdata.sh` [Descarga de datos de RNASeq y miRNASeq; Gráfica de miRNAs y de Archivos]

## Control de calidad y normalizacion (se usan los scripts de R)
   1. `bash sh/DataProcess.sh` [Preprocesamiento de archivos, Pre-Control de Calidad, Plots, Normalización, Post-Control de Calidad]  



# Árbol de Directorios

   - Data
        - Adeno
           - NAD
              - miRNA
                 - 0ce7bbc4-8258-4e88-afbd-26d99b8af1e1.mirbase21.mirnas.quantification.txt (_file_)
                 - ...
              - RNASeq
                 - 1aeb8b9e-bf79-415d-a8dd-a63c3fbe2bab.htseq.counts (_file_)
                 - ...
              - NAD-cases.tsv (_file_)
              - NAD-mirna.tsv (_file_)
           - TAD
              - miRNA
                 - 0a5aa815-9305-4ed3-a0fb-7b6dc6fc503d.mirbase21.mirnas.quantification.txt (_file_)
                 - ...
              - RNASeq
                 - 0b07a88f-e7c7-4fcc-8b95-6a7674c436cc.htseq.counts (_file_)
                 - ...
              - TAD-cases.tsv (_file_)
              - TAD-mirna.tsv (_file_)
           - rdata
              - plots
                 - QC_PRE
                    - 01-biodetection.Rd_01.png (_file_)
                    - 01-biodetection.Rd_02.png (_file_)
                    - 02-countsbio.png (_file_)
                    - 02-protein_coding_barplot.png (_file_)
                    - 02-protein_coding_boxplot.png (_file_)
                    - 03-protein_coding_boxplot_group.png (_file_)
                    - 04-protein_coding_barplot_group.png (_file_)
                    - 05-Lengthbias.png (_file_)
                    - 06-GCbias.png (_file_)
                    - 08-PCAVariance_raw.pdf (_file_)
                    - 09-PCALoading_raw.pdf (_file_)
                    - 10-PCAScore_raw.pdf (_file_)

              - annot.RData (_file_)
              - CancerRaw.RData (_file_)
              - Length.full_GC.full_Between.tmm_Norm_cpm10_arsyn_all.tsv (_file_)
              - Length.full_GC.full_Between.tmm_Norm_cpm10_arsyn.RData (_file_)
              - Length.full_GC.full_Between.tmm_Norm_cpm10_genelist.txt (_file_)
              - Length.full_GC.full_Between.tmm_Norm_cpm10.RData (_file_)
              - Length.full_GC.full_Between.tmm_Norm.RData (_file_)
              - Mean10_ProteinCoding.RData (_file_)
              - NormalRaw.RData (_file_)
              - Norm_cpm10_arsyn.RData (_file_)
              - norm-NAD.tsv (_file_)
              - norm-TAD.tsv (_file_)
              - RawFull.RData (_file_)
              - raw.tsv (_file_)
              - Targets.tsv (_file_)
        - Squamous
           - NSC
              - miRNA
                 - 0ab8ae70-56d3-4ebd-9094-71b3165faa61.mirbase21.mirnas.quantification.txt (_file_)
                 - ...
              - RNASeq
                 - 1ea50a97-2adb-4252-aec8-2fe94d49b0a0.htseq.counts (_file_)
                 - ...
              - NSC-cases.tsv (_file_)
              - NSC-mirna.tsv (_file_)
           - TSC
              - miRNA
                 - 0ab8ae70-56d3-4ebd-9094-71b3165faa61.mirbase21.mirnas.quantification.txt (_file_)
                 - ...
              - RNASeq
                 - 0aa8097b-26b0-4b82-bfab-a5c52d39cf5e.htseq.counts (_file_)
                 - ...
              - TSC-cases.tsv (_file_)
              - TSC-mirna.tsv (_file_)
           - rdata
              - plots
                 - QC_PRE
                    - 01-biodetection.Rd_01.png (_file_)
                    - 01-biodetection.Rd_02.png (_file_)
                    - 02-countsbio.png (_file_)
                    - 02-protein_coding_barplot.png (_file_)
                    - 02-protein_coding_boxplot.png (_file_)
                    - 03-protein_coding_boxplot_group.png (_file_)
                    - 04-protein_coding_barplot_group.png (_file_)
                    - 05-Lengthbias.png (_file_)
                    - 06-GCbias.png (_file_)
                    - 08-PCAVariance_raw.pdf (_file_)
                    - 09-PCALoading_raw.pdf (_file_)
                    - 10-PCAScore_raw.pdf (_file_)

              - annot.RData (_file_)
              - CancerRaw.RData (_file_)
              - Length.full_GC.full_Between.tmm_Norm_cpm10_arsyn_all.tsv (_file_)
              - Length.full_GC.full_Between.tmm_Norm_cpm10_arsyn.RData (_file_)
              - Length.full_GC.full_Between.tmm_Norm_cpm10_genelist.txt (_file_)
              - Length.full_GC.full_Between.tmm_Norm_cpm10.RData (_file_)
              - Length.full_GC.full_Between.tmm_Norm.RData (_file_)
              - Mean10_ProteinCoding.RData (_file_)
              - NormalRaw.RData (_file_)
              - Norm_cpm10_arsyn.RData (_file_)
              - norm-NAD.tsv (_file_)
              - norm-TAD.tsv (_file_)
              - RawFull.RData (_file_)
              - raw.tsv (_file_)
              - Targets.tsv (_file_)

   - Database
      - NAD.txt (_file_)
      - NSC.txt (_file_)
      - TAD.txt (_file_)
      - TSC.txt (_file_)
   - Graphs
      - no_files.png (_file_)
      - no_mirnas_ALL.png (_file_)
   - json
      - qbyfileid.json (_file_)
      - qbyMIRNA.json (_file_)
   - pipeline
      - biomart-20181212.txt (_file_)
   - py
      - casemirna.py (_file_)
      - file_number.py (_file_)
      - getmirnasCounts.py (_file_)
      - Util.py (_file_)
      - downdata.py (_file_)
      - getcases.py (_file_)
   - R
      - 01-Rectify_rnaSeq.R (_file_)
      - 02-PRE-QC.R (_file_)
      - 03-NORM.R (_file_)
      - 04-POST-QC.R (_file_)
   - sh
      - clean.sh (_file_)
      - DataProcess.sh (_file_)
      - downdata.sh (_file_)
   - README.md (_file_)
