###############################################################################
##Get the Work and Data dir
###############################################################################
setwd("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/")
Sys.umask("003")                                                             #?
DATADIR <- 'pipeline/data/'
RDATA <- paste(DATADIR, "rdata", sep="")
dir.create(RDATA)
cat('Data directory: ', DATADIR, '\n')

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
normdir <- '/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/NAD/miRNA'
cat('Checking normal samples \n')
cases <- read.table("NAD-mirna.tsv", header=T,sep="\t")
#cases$V2 <- gsub(".{3}$", "", cases$V4)
cases$mirna_fname <- paste(normdir,cases$mirna_fname, sep="/")
cat('Normal directory: ', normdir, '\n')
#hola <- read.delim(cases$mirna_fname[1], header=T,sep="\t")
hola2 <- read.delim("/home/kevinml/Documentos/INMEGEN/lung-mirnas/Database/NAD/miRNA/457eb90d-9827-4986-bcc1-053d0a4d43b1.mirbase21.mirnas.quantification.txt", header=T,sep="\t")
hola2 <- hola2[,c(1,2)]
#lol <- list(cases$mirna_fname)
#normal<-sapply(cases$mirna_fname, read.delim)
#normal<-bplapply("/home/kevinml/Documentos/INMEGEN/Lung_miRNAs_kevin/mirnas_NAD/", read.delim, sep="\t", header=T)
normal<-bplapply(cases$mirna_fname, read.delim, sep="\t", header=T)
#, col.names=c("EnsemblID", "raw_counts")


##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(hola2, dim)))
stopifnot(nrow(size)==1)
cat('Normal samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(hola2, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))
cat('Genes in normal samples match positions \n')

##Let's keep only the raw counts
normal<-do.call(cbind, lapply(normal, function(x)x[,"raw_counts"]))
targets<-data.frame(File=cases$V2, ID=paste("Ctr", 1:length(cases$V2), sep=""), CASEID=cases$V1)
colnames(normal)<-targets$ID

##Let's change the annotation 
genes<-do.call(rbind,sapply(genes[,1], strsplit, split=".", fixed=TRUE))
colnames(genes)<-c("EnsemblID", "version")

##Save clean data
normal<-list(Counts=normal, Annot=genes, targets=targets)
save(normal, file=paste(RDATA, "NormalRaw.RData", sep="/"), compress="xz")
cat('NormalRaw.RData saved \n')

###############################################################################
## 1) Read TUMOR Data
##      -Exploring the first sample
##      -Check if all samples have the same size
##      -Check if the genes match positions
##      -Let's keep only the raw counts
##      -Let's change the annotation 
##      -Save clean data
##############################################################################
casedir <- 'Data'
##########################
## TAD
ruta <- paste(casedir, "Basal/RNAseq", sep="/")
cases <- read.table(paste(ruta, "db.csv", sep="/"),header=F,sep=",")
cases$V2 <- gsub(".{3}$", "", cases$V2)
cases$V2 <- paste(ruta,cases$V2, sep="/")
cat('Basal directory: ', ruta, '\n')
BasalData<-bplapply(cases$V2, read.delim, sep="\t", header=F, col.names=c("EnsemblID", "raw_counts"))

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(BasalData, dim)))
stopifnot(nrow(size)==1)
cat('Basal samples have the same size \n')

##Check if the genes match positions
genes<-do.call(cbind,lapply(BasalData, function(x)as.character(x[,1])))
genes<-t(unique(t(genes)))
stopifnot(dim(genes)==c(size[1,1], 1))
cat('Genes in Basal samples match positions \n')

##Let's keep only the raw counts
Basal<-do.call(cbind, lapply(BasalData, function(x)x[,"raw_counts"]))
targets<-data.frame(File=cases$V2, ID=paste("tBa", 1:length(cases$V2), sep=""), CASEID=cases$V1)
colnames(Basal)<-targets$ID

##Let's change the annotation 
genes<-do.call(rbind,sapply(genes[,1], strsplit, split=".", fixed=TRUE))
colnames(genes)<-c("EnsemblID", "version")

##Save clean data
Basal<-list(Counts=Basal, Annot=genes, targets=targets)
save(Basal, file=paste(RDATA, "BasalRaw.RData", sep="/"), compress="xz")
cat('BasalRaw.RData saved \n')

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

names(annot)<-c("EnsemblID", "Chr", "Start", "End", "GC", "Type")
annot$Length<-abs(annot$End - annot$Start)

## Exploring the header
head(annot)

##Remove not conventional chromosomes
annot<-as.data.frame(annot)
levels(annot$Chr)
annot<-annot[annot$Chr%in%c(as.character(1:22), "X", "Y"),]
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
load(file=paste(RDATA, "BasalRaw.RData", sep="/"))
load(file=paste(RDATA, "LumARaw.RData", sep="/"))
load(file=paste(RDATA, "LumBRaw.RData", sep="/"))
load(file=paste(RDATA, "Her2Raw.RData", sep="/"))

##M=normal|tumor
M<-cbind(normal$Counts, Basal$Counts, LumA$Counts, LumB$Counts, Her2$Counts)
cat('Total number of features and samples: ', paste(dim(M), collapse=" ,"), '\n')

##targets=normal+tumor
targets<-rbind(normal$targets, Basal$targets, LumA$targets, LumB$targets, Her2$targets)
targets<-data.frame(targets, stringsAsFactors=FALSE)
targets$File<-as.character(targets$File)
targets$ID<-as.character(targets$ID)
dim(targets)

##check M y targets integrity
stopifnot(nrow(targets)==ncol(M))
cat('Number of counts columns match sample number\n')

## Annot=genes|annot
## check if genes from normal and tumor match
stopifnot(normal$Annot==Basal$Annot)
stopifnot(normal$Annot==LumA$Annot)
stopifnot(normal$Annot==LumB$Annot)
stopifnot(normal$Annot==Her2$Annot)
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

