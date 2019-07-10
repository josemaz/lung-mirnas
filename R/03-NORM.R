###############################################################################
##Get the Tissue Type (Adenocarcinoma or Squamous Cells carcinoma)  
###############################################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Execution: Rscript --vanilla R/03-NORM.R (tissue type); tissue type: Adeno or Squamous", call.=FALSE)
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
##Usefull Libraries
###############################################################################
library("NOISeq")
library("BiocParallel")
library("EDASeq")

register(MulticoreParam(workers=detectCores()-1, progress=TRUE))#Linux

###############################################################################
## Filter genes with low expression
cat("#####################\n")
cat("Step 3: Normalization\n")
cat("#####################\n")

RDATA <- paste("Data",Tissue, "rdata", sep = "/")
PLOTSDIR <-paste(RDATA, "plots", sep = "/")
w <- 1024
h <- 1024
p <- 24

load(file=paste(RDATA, "RawFull.RData", sep="/"))
### We keep only genes with mean expression count > 10 
exp.genes <- apply(full$M, 1, function(x) mean(x)>10)
egtable <- table(exp.genes)
cat("There are", egtable[[1]], "genes with mean expression count < 10", egtable[[2]], "with mean count > 10 \n")

##Filtering low expression genes
mean10 <- list(M=full$M[exp.genes, ], Annot=full$Annot[exp.genes, ], Targets=full$Targets)
rownames(mean10$Annot) <- rownames(mean10$M)

cat("Filtering protein coding features \n")
protein.coding <- rownames(mean10$Annot[mean10$Annot$Type == "protein_coding", ])
cat("There are ", length(protein.coding), " features with type: protein_coding \n")
mean10 <- list(M=mean10$M[protein.coding, ], Annot=mean10$Annot[protein.coding, ], Targets=mean10$Targets)

# Escribe los raw counts (filtrados por biotipo y mean > 10)
# esto para usar con otras normalizaciones que no esten en este trabajo
rawM <- cbind(gene=as.character(mean10$Annot$EnsemblID), mean10$M)
write.table(rawM, file=paste(RDATA, "raw.tsv", sep="/"), sep="\t", quote=FALSE, row.names=FALSE)

cat("\nSaving Mean10_ProteinCoding.RData...\n") 
save(mean10, file=paste(RDATA, "Mean10_ProteinCoding.RData", sep="/"), compress="xz")
cat("Saved.\n")

##########################################
## Normalization: GC, Lenght, Tmm
##########################################
load(file=paste(RDATA, "Mean10_ProteinCoding.RData", sep="/"))
ln.data <- withinLaneNormalization(mean10$M, mean10$Annot$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data, mean10$Annot$GC, which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
step1 <- "Length.full"
step2 <- "GC.full"
step3 <- "Between.tmm"

##########################################
## ARSyN to reduce batch effect
##########################################
norm.data <- mean10
norm.data$M <- norm.counts

cat("\nSaving", paste(step1, step2, step3, "Norm.RData ...", sep = "_"), "\n")
save(norm.data, file=paste(RDATA, paste(step1, step2, step3, "Norm.RData", sep = "_"), sep="/"), compress="xz")
cat("Saved.\n")

## LOW EXPRESSION FILTERING
### ¿Vuelve a filtrar?
norm.data.cpm10 <- filtered.data(norm.data$M, factor=norm.data$Targets$Group, 
                                 norm=TRUE, cv.cutoff=100, cpm=10)
filtered <- nrow(norm.data$M)-nrow(norm.data.cpm10)
cat("\nAfter normalization. There are", filtered, "genes with counts per million mean < 10", 
    nrow(norm.data.cpm10), "with counts per million mean > 10 \n")

norm.data.cpm10 <-list(M = norm.data.cpm10, 
                       Annot = norm.data$Annot[row.names(norm.data$Annot) %in% row.names(norm.data.cpm10),],
                       Targets = norm.data$Targets)

stopifnot(nrow(norm.data.cpm10$M) == nrow(norm.data.cpm10$Annot))
stopifnot(row.names(norm.data.cpm10$M) == row.names(norm.data.cpm10$Annot))

cat("\nSaving", paste(step1, step2, step3, "Norm_cpm10.RData ...", sep = "_"), "\n")
save(norm.data.cpm10, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10.RData", sep = "_"), sep="/"), compress="xz")
cat("Saved.\n")

mydata <- NOISeq::readData(
  data = norm.data.cpm10$M, 
  factors = norm.data.cpm10$Targets[, "Group",drop = FALSE])

cat("\nPerforming ARSyN for batch correction\n")
myARSyN <- ARSyNseq(mydata, norm = "n", logtransf = FALSE)

cat('\nARSyN data. Final dimensions: ', paste(dim(assayData(myARSyN)$exprs), collapse=", "), '.\n')

##Saving everything
norm.data_cpm10_arsyn <- list(M = assayData(myARSyN)$exprs, Annot = norm.data.cpm10$Annot, 
                              Targets = norm.data.cpm10$Targets)

stopifnot(nrow(norm.data_cpm10_arsyn$M) == nrow(norm.data_cpm10_arsyn$Annot))
stopifnot(all(row.names(norm.data_cpm10_arsyn$M) == row.names(norm.data_cpm10_arsyn$Annot)))

cat("\nSaving", paste(step1, step2, step3, "Norm_cpm10_arsyn.RData ...", sep = "_"), "\n")
save(norm.data_cpm10_arsyn, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_arsyn.RData", sep = "_"), sep="/"), compress="xz")
cat("Saved.")

cat("\nGenerating data matrices with arsyn for Aracne...\n")
## Data matrices for Aracne
## ALL = healthy | cancer
M <- as.data.frame(norm.data_cpm10_arsyn$M)
M <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), M)

#healthy samples
NAD_NSC <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == normalTissue])
NAD_NSC <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), NAD_NSC)
#tumor samples
TAD_TSC <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == cancerTissue])
TAD_TSC <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), TAD_TSC)

#EnsemblIDs
symbols <-as.character(norm.data_cpm10_arsyn$Annot$EnsemblID)

cat("\nSaving arsyn data\n")

cat("Saving", paste(step1, step2, step3, "Norm_cpm10_arsyn_all.tsv", sep = "_"), "\n")
write.table(M, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_arsyn_all.tsv", sep = "_"), sep="/"), 
            sep="\t", quote=FALSE, row.names=FALSE)
cat("Saving", paste(step1, step2, step3, "Norm_cpm10_genelist.txt", sep = "_"), "\n")
write.table(symbols, file = paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_genelist.txt", sep = "_"), sep="/"), 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("arsyn data saved!\n")

cat("\nSaving norm-", normalTissue, ".tsv", "\n")
write.table(NAD_NSC, file =paste(RDATA, "/norm-", normalTissue, ".tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
cat("File saved.\n")
cat("\nSaving norm-", cancerTissue, ".tsv", "\n")
write.table(TAD_TSC, file =paste(RDATA, "/norm-", cancerTissue, ".tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
cat("File saved.\n")

cat("\nSaving", "Norm_cpm10_arsyn.RData", "\n")
save(norm.data_cpm10_arsyn, file=paste(RDATA,"Norm_cpm10_arsyn.RData", sep="/"), compress="xz")
cat("File saved.\n")

cat("\nEnd of normalization texting\n")


