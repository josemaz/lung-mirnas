###############################################################################
##Get the Tissue Type (Adenocarcinoma or Squamous Cells carcinoma)  
###############################################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Execution: Rscript --vanilla R/miRNA-Normlize.R (tissue type); tissue type: Adeno or Squamous", call.=FALSE)
} 
Tissue <- args[1] 
if (Tissue == "Adeno"){
  normalTissue <- "NAD"
  cancerTissue <- "TAD"
} else if (Tissue == "Squamous"){
  normalTissue <- "NSC"
  cancerTissue <- "TSC"  
}
###############################################################################
##Get the Work and Data dir
###############################################################################
Sys.umask("003")
RDATA <- paste("Data",Tissue, "rdata", "miRNA", sep = "/")
dir.create(RDATA)
cat('Data directory: ', "Data/", '\n')
###############################################################################
##Usefull Libraries
###############################################################################
library("BiocParallel")
library("NOISeq")
###################################################
cat("######################\n")
cat(" Step 1: Rectify Data \n")
cat("######################\n")
casedir <- paste("Data",Tissue, sep = "/")
########################## 
## Normal NAD OR NSC
ruta <- paste(casedir, normalTissue, "miRNA", sep="/")
cases <- read.table(paste(ruta,"index.txt", sep ="/"), header=F, sep='\t')

cases$V2 <- paste(ruta,cases$V2, sep="/")
cat('Normal directory: ', ruta, '\n')
NormalData <-bplapply(cases$V2, read.delim, sep="\t", header=T,
                     col.names=c("miRNA_ID", "read_counts","RPM_mapped","cross"))

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(NormalData, dim)))
stopifnot(nrow(size)==1)
cat('Normal samples have the same size \n')

##Get and Check if the miRNAs ID with samples
miRNAs<-do.call(cbind,lapply(NormalData, function(x)as.character(x[,1])))
miRNAs<-t(unique(t(miRNAs)))
stopifnot(dim(miRNAs)==c(size[1,1], 1))
cat('miRNAs in Normal samples are OK \n')

##Lets keep only the raw counts
rcounts<-do.call(cbind, lapply(NormalData, function(x)x[,"read_counts"]))
targets<-data.frame(File=cases$V2, ID=paste(normalTissue, 1:length(cases$V2), sep=""), CASEID=cases$V1)
colnames(rcounts)<-targets$ID

Normal<-list(Counts=rcounts, Annot=miRNAs, targets=targets)
save(Normal, file=paste(RDATA, "NormalRaw.RData", sep="/"), compress="xz")
cat('NormalRaw.RData saved \n')


##########################
## Cancer TAD OR TSC
ruta <- paste(casedir, cancerTissue, "miRNA", sep="/")
cases <- read.table(paste(ruta,"index.txt", sep ="/"), header=F, sep='\t')
cases$V2 <- paste(ruta,cases$V2, sep="/")
cat('Cancer directory: ', ruta, '\n')
CancerData <-bplapply(cases$V2, read.delim, sep="\t", header=T,
                     col.names=c("miRNA_ID", "read_counts","RPM_mapped","cross"))

##Check if all samples have the same size
size<-unique(do.call(rbind,lapply(CancerData, dim)))
stopifnot(nrow(size)==1)
cat('Cancer samples have the same size \n')

##Get and Check if the miRNAs ID with samples
miRNAs<-do.call(cbind,lapply(CancerData, function(x)as.character(x[,1])))
miRNAs<-t(unique(t(miRNAs)))
stopifnot(dim(miRNAs)==c(size[1,1], 1))
cat('miRNAs in Cancer samples are OK \n')

##Lets keep only the raw counts
rcounts<-do.call(cbind, lapply(CancerData, function(x)x[,"read_counts"]))
targets<-data.frame(File=cases$V2, ID=paste(cancerTissue, 1:length(cases$V2), sep=""), CASEID=cases$V1)
colnames(rcounts)<-targets$ID

Cancer<-list(Counts=rcounts, Annot=miRNAs, targets=targets)
save(Cancer, file=paste(RDATA, "CancerRaw.RData", sep="/"), compress="xz")
cat('CancerRaw.RData saved \n')


##############################################################################
## Merging count and annotation
cat('Merging counts and annotations \n')
load(file=paste(RDATA, "CancerRaw.RData", sep="/"))
load(file=paste(RDATA, "NormalRaw.RData", sep="/"))

##M = Normal|Cancer
M<-cbind(Normal$Counts, Cancer$Counts)
# M<-cbind(Normal$Counts, LumA$Counts)
cat('Total number of features and samples: ', paste(dim(M), collapse=" ,"), '\n')

##targets=normal+tumor
targets<-rbind(Normal$targets, Cancer$targets)
# targets<-rbind(Normal$targets, LumA$targets)
targets<-data.frame(targets, stringsAsFactors=FALSE)
targets$File<-as.character(targets$File)
targets$ID<-as.character(targets$ID)
dim(targets)

##check M y targets integrity
stopifnot(nrow(targets)==ncol(M))
cat('Number of counts columns match sample number\n')

## Annot=miRNAs|annot
## check if miRNAs from normal and tumor match
stopifnot(Normal$Annot == Cancer$Annot)
cat('miRNA IDs from samples match \n')

Annot<-Normal$Annot
Annot<-data.frame(Annot, stringsAsFactors=FALSE)
Annot$Row<-1:nrow(Annot)
head(Annot)

##Save the clean Data
full<-list(M=M, Annot=Annot, Targets=targets)
cat("Saving RawFull.RData...\n")
save(full, file=paste(RDATA, "RawFull.RData", sep="/"), compress="xz")
cat("RawFull.RData saved.\n")


###################################################
cat("\n#################\n")
cat("   Step 2: QC\n")
cat("#################\n")

load(file=paste(RDATA, "RawFull.RData", sep="/"))
full$Targets$Group<-factor(substr(full$Targets$ID, start=1, stop=3))
f <- full$Targets[, "Group",drop=FALSE]
rownames(f) <- full$Targets$ID

# establecer el umbral para que si hay menos de 382 ceros 
# (~ mitad del total de muestras) sobre el mirna, entonces borrarlo
# Validado por el Dr. Jesus Espinal.
mitad <- apply(full$M == 0, 1, sum) < 382
full$M <- full$M[mitad, ]
full$Annot <- full$Annot[mitad,]
full$Annot$Row <- 1:nrow(full$Annot)

mydata <- NOISeq::readData(data = full$M, factors = f)
myPCA = dat(mydata, type = "PCA")
par(mfrow = c(1,2))
explo.plot(myPCA, factor = "Group")

myARSyN = ARSyNseq(mydata, factor = "Group", batch = FALSE, norm = "rpkm",  logtransf = FALSE)
myPCA = dat(myARSyN, type = "PCA")
explo.plot(myPCA, factor = "Group")
# heatmap(full$M)

cat('ARSyN data. Final dimensions: ', paste(dim(assayData(myARSyN)$exprs), collapse=", "), '.\n')

data_cpm10_arsyn <- list(M = assayData(myARSyN)$exprs, Annot = full$Annot, 
                         Targets = full$Targets)
rownames(data_cpm10_arsyn$M) <- row.names(data_cpm10_arsyn$Annot)
stopifnot(nrow(data_cpm10_arsyn$M) == nrow(data_cpm10_arsyn$Annot))

M <- as.data.frame(data_cpm10_arsyn$M)
M <- cbind(mirna=as.character(data_cpm10_arsyn$Annot$Annot), M)

write.table(M, file = paste(RDATA, "all.tsv", sep="/"),
            sep="\t", quote=FALSE, row.names=FALSE)

normal <- as.data.frame(data_cpm10_arsyn$M[,data_cpm10_arsyn$Targets$Group == normalTissue])
normal <- cbind(mirna=as.character(data_cpm10_arsyn$Annot$Annot), normal)
write.table(normal, file = paste(RDATA, "Normal.tsv", sep="/"), 
            sep="\t", quote=FALSE, row.names=FALSE)

cancer <- as.data.frame(data_cpm10_arsyn$M[,data_cpm10_arsyn$Targets$Group == cancerTissue])
cancer <- cbind(mirna=as.character(data_cpm10_arsyn$Annot$Annot), cancer)
write.table(cancer, file = paste(RDATA, "Cancer.tsv", sep="/"),
            sep="\t", quote=FALSE, row.names=FALSE)

all <-
  
  Targets <- as.data.frame(data_cpm10_arsyn$Targets)
write.table(Targets, file = paste(RDATA, "Targets.tsv", sep="/"),
            sep="\t", quote=FALSE, row.names=FALSE)

