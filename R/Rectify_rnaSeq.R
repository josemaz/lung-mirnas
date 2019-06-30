###############################################################################
##Get the Work and Data dir
###############################################################################
#ruta donde este la carpeta Data
Sys.umask("003")
DATADIR <- 'pipeline/data/'
RDATA <- paste(DATADIR, "rdata", sep="")
dir.create(RDATA)
cat('Data directory: ', DATADIR, '\n') #concatena e imprime

###############################################################################
##Usefull Libraries
###############################################################################
library("BiocParallel")
library("parallel")
#register(SnowParam(workers=detectCores()-1, progress=TRUE))#Windows
register(MulticoreParam(workers=detectCores()-1, progress=TRUE))#Linux

###############################################################################
## 1) Read Normal Data
##      -Exploring the first sample
##      -Check if all samples have the same size
##      -Check if the genes match positions
##      -Let's keep only the raw counts
##      -Let's change the annotation 
##      -Save clean data
##############################################################################

##########################
## Control - NAD
normdir <- 'Data/NAD/RNAseq/' #ruta donde este el archivo index.txt
cat('Checking normal samples \n')
cases <- read.table(paste(normdir,"index.txt", sep =""), header=F, sep='\t') #index.txt ahora es cases
cases$V2 <- gsub(".{3}$", "", cases$V2) #Sustitucion de todos los matches de un string
cases$V2 <- paste(normdir,cases$V2, sep="")
cat('Normal directory: ', normdir, '\n')
normal<-bplapply(cases$V2 , read.delim, sep="\t", header=F, col.names=c("EnsemblID", "raw_counts")) #bplapply - Aplica la funcion read.delim para cada elemento de cases$ V2
#read-delim - lee archivos de texto donde filas-casos & columnas-muesstras

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(normal, dim)))
stopifnot(nrow(size)==1)
cat('Normal samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(normal, function(x)as.character(x[,1]))) #genes es una matriz con todos los ensemble ID
genes<-t(unique(t(genes))) #Se camnia el formato de la matriz, ya es mas una lista que una matriz
stopifnot(dim(genes)==c(size[1,1], 1)) #Checar que se tengan ciertas dimensiones 
cat('Genes in normal samples match positions \n')

##Let's keep only the raw counts
normal<-do.call(cbind, lapply(normal, function(x)x[,"raw_counts"])) #Como los ensemble Id ya estan en genes, solo nos quedamos con las cuentas en "normal"
targets<-data.frame(File=cases$V2, ID=paste("Ctr", 1:length(cases$V2), sep=""), CASEID=cases$V1) #File = el path de los archivos
#ID = numero de muestra/archivo, CaseID = CaseID
colnames(normal)<-targets$ID #normal tiene muchas columnas, donde cada columna es un gen y sus filas son las cuentas de un transcrito

##Let's change the annotation 
genes<-do.call(rbind,sapply(genes[,1], strsplit, split=".", fixed=TRUE))
colnames(genes)<-c("EnsemblID", "version") # genes es un dataframe de 2 columnas, en la primera esta el ensemble ID y en la otra la verison
write.table(genes[,"EnsemblID"], "ensembleIDs_Normal.tsv", sep="\t", row.names = FALSE) #Se hace un archivo con los ensemble IDs

##Save clean data
normal<-list(Counts=normal, Annot=genes, targets=targets)
save(normal, file=paste(RDATA, "NormalRaw.RData", sep="/"), compress="xz")
cat('NormalRaw.RData saved \n')
#Al final normal es una lista donde en el primer elemento estan las cuentas de cada gen "Ctr", en el segundo elemento estan los ensembleID, en el tercer elemento estan
#en el tercer elemento estan los paths de los archivos, en el cuarto los "Ctr" y en el quinto los CaseID
###############################################################################
## 1) Read TUMOR Data
##      -Exploring the first sample
##      -Check if all samples have the same size
##      -Check if the genes match positions
##      -Let's keep only the raw counts
##      -Let's change the annotation 
##      -Save clean data
##############################################################################
casedir <- 'Data/TAD/RNAseq/' #ruta donde este el archivo index.txt
##########################
## TAD
cases <- read.table(paste(casedir, "index.txt", sep=""),header=F,sep="\t")
cases$V2 <- gsub(".{3}$", "", cases$V2)
cases$V2 <- paste(casedir,cases$V2, sep="")
cat('Cancerous directory: ', casedir, '\n')
CancerData<-bplapply(cases$V2, read.delim, sep="\t", header=F, col.names=c("EnsemblID", "raw_counts"))

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(CancerData, dim)))
stopifnot(nrow(size)==1)
cat('Cancer samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(CancerData, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))
cat('Genes in cancer samples match positions \n')

##Let's keep only the raw counts
cancer<-do.call(cbind, lapply(CancerData, function(x)x[,"raw_counts"]))
targets<-data.frame(File=cases$V2, ID=paste("tBa", 1:length(cases$V2), sep=""), CASEID=cases$V1)
colnames(cancer)<-targets$ID

##Let's change the annotation 
genes<-do.call(rbind,sapply(genes[,1], strsplit, split=".", fixed=TRUE))
colnames(genes)<-c("EnsemblID", "version")
write.table(genes[,"EnsemblID"], "ensembleIDs_Cancer.tsv", sep="\t", row.names = FALSE)

##Save clean data
cancer<-list(Counts=cancer, Annot=genes, targets=targets)
save(cancer, file=paste(RDATA, "CancerRaw.RData", sep="/"), compress="xz")
cat('CancerRaw.RData saved \n')

##############################################################################
## 3) Read Biomart Data
##      -Read the file
##      -Exploring the header
##      -Remove not conventional chromosomes
##      -Keep only annotated EntrezID individuals
##      -Remove the ones without Symbol | HGNCID
##      -Save clean data
###############################################################################
## Read the file
cat('Working with annotation file: biomart-20181212.txt \n')
annot<-read.delim(file="pipeline/biomart-20181212.txt", sep="\t")

#names(annot)<-c("EnsemblID", "Start", "End", "GC", "Type", "Chr", "Band")
names(annot)<-c("EnsemblID", "Chr", "Start", "End", "GC", "Type")
annot$Length<-abs(annot$End - annot$Start)

## Exploring the header
head(annot)

##Remove not conventional chromosomes
annot<-as.data.frame(annot)
levels(annot$Chr) #arroja los valores unicos, cromosomas unicos
annot<-annot[annot$Chr%in%c(as.character(1:22), "X", "Y"),] #Solo se dejan aquellos cromosomas del1 al 22 y los sexuales
annot$Chr<-droplevels(annot$Chr)
cat('Non conventional chromosomes removed \n')

uniq.annot <- length(unique(annot$EnsemblID)) == nrow(annot)
if(uniq.annot) {
  cat('Unique EnsemblIDs in annotation file\n')
} else {
  cat('Repeated EnsemblIDs in annotation file. \n')
  stop()
}

cat('Annotation file. Final dimension: ', paste(dim(annot), collapse=", "), '\n')
## Save clean data
save(annot, file=paste(RDATA, "annot.RData", sep="/"), compress="xz")
cat('annot.RData saved \n')

##############################################################################
##4) Merging count and annotation
##	    -M=normal|tumour
##      -targets=normal+tumor
##          -Check M y targets integrity
##	    -Annot=genes|annot
##          -Check if genes from normal and tumour match
##          -Are there repeated EntrezID?
##          -Add Biomart data
##          -Are there duplicated IDs?
##              -Keep the matching Symbol.x==Symbol.y
##              -Keep the matching Symbol.x==Symbol.y
##              -Keep the one with lowest GC content
##          -Did we lost any row in the process?
##          -Save the clean Data
##############################################################################
cat('Merging counts and annotations \n')
load(file=paste(RDATA, "annot.RData", sep="/"))
load(file=paste(RDATA, "NormalRaw.RData", sep="/"))
load(file=paste(RDATA, "CancerRaw.RData", sep="/"))

##M=normal|tumor
M<-cbind(normal$Counts, cancer$Counts)
cat('Total number of features and samples: ', paste(dim(M), collapse=" ,"), '\n')

##targets=normal+tumor
targets<-rbind(normal$targets, cancer$targets)
targets<-data.frame(targets, stringsAsFactors=FALSE)
targets$File<-as.character(targets$File)
targets$ID<-as.character(targets$ID)
dim(targets)

##check M y targets integrity
stopifnot(nrow(targets)==ncol(M))
cat('Number of counts columns match sample number\n')

## Annot=genes|annot
## check if genes from normal and tumor match
stopifnot(normal$Annot==cancer$Annot)

cat('Genes from normal and tumor samples match \n')

Annot<-normal$Annot
##Are there repeated EnsemblID?
stopifnot(length(unique(Annot[, "EnsemblID"]))==nrow(Annot))
cat('No duplicated EnsemblIDs\n')

Annot<-data.frame(Annot, stringsAsFactors=FALSE)
Annot$Row<-1:nrow(Annot) ##Just to maintain the original order
head(Annot)

cat('Adding biomart data\n')
Annot<-merge(x=Annot, y=annot, by.x="EnsemblID", by.y="EnsemblID",
             all=FALSE, all.x=TRUE, all.y=FALSE, sort=FALSE)
cat('Merged file. Final dimensions: ', paste(dim(Annot), collapse=", "), '.\n')
dim(M)
extra.rows <- nrow(Annot)-nrow(M) 
cat('There are ', extra.rows, ' extra rows in the counts matrix.\n')

##Are there duplicated IDs?
Annot<-Annot[order(Annot$Row, Annot$Chr),]
ids<-head(which(duplicated(Annot$Row)))
ids
if (identical(ids, integer(0))) {
  cat('There are no duplicated IDs \n')
} else {
  cat('Si hay duplicados copiar codigo del pipeline')
}

############### QUE PASA EN ESTOS CASOS? #########
gctable <- table(is.na(Annot$GC))
cat("There are", gctable[[1]], "entries with GC info and", gctable[[2]], "GC entries missing \n")

##Did we lost any row in the process?
stopifnot(all(Annot$Row %in% 1:nrow(M)))
stopifnot(all(1:nrow(M) %in% Annot$Row))

##Save the clean Data
full<-list(M=M, Annot=Annot, Targets=targets)
save(full, file=paste(RDATA, "RawFull.RData", sep="/"), compress="xz")
cat("Saving RawFull.RData \n")

##############################################################################
## GREAT JOB!!! YOU MADE IT TILL THE END!!!!
###############################################################################
cat("###########################\n")
cat("End of Step 1: Get The Data\n")
cat("###########################\n")

