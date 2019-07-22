###############################################################################
##Create all Expression Matrices
###############################################################################

# Usefull libraries
library("BiocParallel")
#################################
# Se lee la informacion de Adeno
#################################
#Se carga el archivo donde viene la relacion muestra (NAD1, NAD2, ...) con CaseID (TCGA-38-4625, ...)
#Leemos los RNAs
load("Data/Adeno/rdata/RNA/RawFull.RData")
caseid_RNA_Adeno <- full$Targets$CASEID
muestra_RNA_Adeno <- full$Targets$ID

#VAMOS A LEER LOS miRNAs
load("Data/Adeno/rdata/miRNA/RawFull.RData")
caseid_miRNA_Adeno <- full$Targets$CASEID
muestra_miRNA_Adeno <- full$Targets$ID

tablita <- data.frame(caseid_RNA_Adeno,caseid_miRNA_Adeno, muestra_RNA_Adeno, muestra_miRNA_Adeno)
#Si todos los case_id son iguales (verdaderos) se añaden las muestras abajo
if ( all(caseid_RNA_Adeno == caseid_miRNA_Adeno, na.rm = FALSE) & all(muestra_RNA_Adeno == muestra_miRNA_Adeno, na.rm = FALSE) ){
  #Pegamos las tablitas de NAD
  cat("\nGenerating Expresion matrix of NAD...")
  RNA_NAD <- read.table("Data/Adeno/rdata/RNA/norm-NAD.tsv", sep="\t", header=T )
  miRNA_NAD <- read.table("Data/Adeno/rdata/miRNA/Normal.tsv", sep="\t", header=T )
  EM_NAD_full <- rbind(RNA_NAD[,-1],miRNA_NAD[,-1])
  row.names(EM_NAD_full) <- c(as.character(RNA_NAD[,1]) , as.character(miRNA_NAD[,1]))
  write.table(EM_NAD_full, file = "Data/Adeno/NAD/EM_RNA_miRNA_NAD.tsv", sep = "\t")
  cat("\nFile saved in: Data/Adeno/NAD/")
  cat("\nDone.\n")
  #Pegamos las tablitas de TAD
  cat("\nGenerating Expresion matrix of TAD...")
  RNA_TAD <- read.table("Data/Adeno/rdata/RNA/norm-TAD.tsv", sep="\t", header=T )
  miRNA_TAD <- read.table("Data/Adeno/rdata/miRNA/Cancer.tsv", sep="\t", header=T)
  EM_TAD_full <- rbind(RNA_TAD[,-1],miRNA_TAD[,-1])  
  row.names(EM_TAD_full) <- c(as.character(RNA_TAD[,1]) , as.character(miRNA_TAD[,1]))
  write.table(EM_TAD_full, file = "Data/Adeno/TAD/EM_RNA_miRNA_TAD.tsv", sep = "\t")
  cat("\nFile saved in: Data/Adeno/TAD/")
  cat("\nDone.\n")
}

###################################
# Se lee la informacion de Squamous
###################################
#Se carga el archivo donde viene la relacion muestra (NAD1, NAD2, ...) con CaseID (TCGA-38-4625, ...)
#Leemos los RNAs
load("Data/Squamous/rdata/RNA/RawFull.RData")
caseid_RNA_Squamous <- full$Targets$CASEID
muestra_RNA_Squamous <- full$Targets$ID

#VAMOS A LEER LOS miRNAs
load("Data/Squamous/rdata/miRNA/RawFull.RData")
caseid_miRNA_Squamous <- full$Targets$CASEID
muestra_miRNA_Squamous <- full$Targets$ID

tablita <- data.frame(caseid_RNA_Squamous,caseid_miRNA_Squamous, muestra_RNA_Squamous, muestra_miRNA_Squamous)
#Si todos los case_id son iguales (verdaderos) se añaden las muestras abajo
if ( all(caseid_RNA_Squamous == caseid_miRNA_Squamous, na.rm = FALSE) & all(muestra_RNA_Squamous == muestra_miRNA_Squamous, na.rm = FALSE) ){
  #Pegamos las tablitas de NSC
  cat("\nGenerating Expresion matrix of NSC...")
  RNA_NSC <- read.table("Data/Squamous/rdata/RNA/norm-NSC.tsv", sep="\t", header=T )
  miRNA_NSC <- read.table("Data/Squamous/rdata/miRNA/Normal.tsv", sep="\t", header=T )
  EM_NSC_full <- rbind(RNA_NSC[,-1],miRNA_NSC[,-1])
  row.names(EM_NSC_full) <- c(as.character(RNA_NSC[,1]) , as.character(miRNA_NSC[,1]))
  write.table(EM_NSC_full, file = "Data/Squamous/NSC/EM_RNA_miRNA_NSC.tsv", sep = "\t")
  cat("\nFile saved in: Data/Squamous/NSC/")
  cat("\nDone.\n")
  #Pegamos las tablitas de TSC
  cat("\nGenerating Expresion matrix of TSC...")
  RNA_TSC <- read.table("Data/Squamous/rdata/RNA/norm-TSC.tsv", sep="\t", header=T )
  miRNA_TSC <- read.table("Data/Squamous/rdata/miRNA/Cancer.tsv", sep="\t", header=T )
  EM_TSC_full <- rbind(RNA_TSC[,-1],miRNA_TSC[,-1])  
  row.names(EM_TSC_full) <- c(as.character(RNA_TSC[,1]) , as.character(miRNA_TSC[,1]))
  write.table(EM_TSC_full, file = "Data/Squamous/TSC/EM_RNA_miRNA_TSC.tsv", sep = "\t")
  cat("\nFile saved in: Data/Squamous/TSC/")
  cat("\nDone.\n\n")
}

