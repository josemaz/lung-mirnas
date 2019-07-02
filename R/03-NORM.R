library("NOISeq")
library("BiocParallel")
library("EDASeq")

register(MulticoreParam(workers=detectCores()-1, progress=TRUE))#Linux

###############################################################################
## Filter genes with low expression
cat("#################\n")
cat("Step 3: Normalization\n")
cat("#################\n")

DATADIR <- 'pipeline/data/'
RDATA <- paste(DATADIR, "rdata", sep = "")
PLOTSDIR <-paste(DATADIR, "plots", sep = "")
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

cat("Saving Mean10_ProteinCoding.RData \n") 
save(mean10, file=paste(RDATA, "Mean10_ProteinCoding.RData", sep="/"), compress="xz")

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

cat("Saving", paste(step1, step2, step3, "Norm.RData", sep = "_"), "\n")
save(norm.data, file=paste(RDATA, paste(step1, step2, step3, "Norm.RData", sep = "_"), sep="/"), compress="xz")

## LOW EXPRESSION FILTERING
### ¿Vuelve a filtrar?
norm.data.cpm10 <- filtered.data(norm.data$M, factor=norm.data$Targets$Group, 
                                 norm=TRUE, cv.cutoff=100, cpm=10)
filtered <- nrow(norm.data$M)-nrow(norm.data.cpm10)
cat("After normalization. There are", filtered, "genes with counts per million mean < 10", 
    nrow(norm.data.cpm10), "with counts per million mean > 10 \n")

norm.data.cpm10 <-list(M = norm.data.cpm10, 
                       Annot = norm.data$Annot[row.names(norm.data$Annot) %in% row.names(norm.data.cpm10),],
                       Targets = norm.data$Targets)

stopifnot(nrow(norm.data.cpm10$M) == nrow(norm.data.cpm10$Annot))
stopifnot(row.names(norm.data.cpm10$M) == row.names(norm.data.cpm10$Annot))

cat("Saving", paste(step1, step2, step3, "Norm_cpm10.RData", sep = "_"), "\n")
save(norm.data.cpm10, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10.RData", sep = "_"), sep="/"), compress="xz")

mydata <- NOISeq::readData(
  data = norm.data.cpm10$M, 
  factors = norm.data.cpm10$Targets[, "Group",drop = FALSE])

cat("Performing ARSyN for batch correction")
myARSyN <- ARSyNseq(mydata, norm = "n", logtransf = FALSE)

cat('ARSyN data. Final dimensions: ', paste(dim(assayData(myARSyN)$exprs), collapse=", "), '.\n')

##Saving everything
norm.data_cpm10_arsyn <- list(M = assayData(myARSyN)$exprs, Annot = norm.data.cpm10$Annot, 
                              Targets = norm.data.cpm10$Targets)

stopifnot(nrow(norm.data_cpm10_arsyn$M) == nrow(norm.data_cpm10_arsyn$Annot))
stopifnot(all(row.names(norm.data_cpm10_arsyn$M) == row.names(norm.data_cpm10_arsyn$Annot)))

cat("Saving", paste(step1, step2, step3, "Norm_cpm10_arsyn.RData", sep = "_"), "\n")
save(norm.data_cpm10_arsyn, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_arsyn.RData", sep = "_"), sep="/"), compress="xz")

cat("Generating data matrices with arsyn for Aracne\n")
## Data matrices for Aracne
## ALL = healthy | cancer
M <- as.data.frame(norm.data_cpm10_arsyn$M)
M <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), M)

#tumor samples
Ctr <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == "Ctr"])
Ctr <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), Ctr)
St1 <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == "s1-"])
St1 <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), St1)
St2 <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == "s2-"])
St2 <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), St2)
St3 <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == "s3-"])
St3 <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), St3)
# St4 <- as.data.frame(norm.data_cpm10_arsyn$M[,norm.data_cpm10_arsyn$Targets$Group == "s4-"])
# St4 <- cbind(gene=as.character(norm.data_cpm10_arsyn$Annot$EnsemblID), St3)

#EnsemblIDs
symbols <-as.character(norm.data_cpm10_arsyn$Annot$EnsemblID)

cat("Saving arsyn data\n")

cat("Saving", paste(step1, step2, step3, "Norm_cpm10_arsyn_all.tsv", sep = "_"), "\n")
write.table(M, file=paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_arsyn_all.tsv", sep = "_"), sep="/"), 
            sep="\t", quote=FALSE, row.names=FALSE)
cat("Saving", paste(step1, step2, step3, "Norm_cpm10_genelist.txt", sep = "_"), "\n")
write.table(symbols, file = paste(RDATA, paste(step1, step2, step3, "Norm_cpm10_genelist.txt", sep = "_"), sep="/"), 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

cat("Saving", "norm-Control.tsv", "\n")
write.table(Ctr, file =paste(RDATA,"norm-Control.tsv", sep="/"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Saving", "norm-Stage1.tsv", "\n")
write.table(St1, file =paste(RDATA,"norm-Stage1.tsv", sep="/"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Saving", "norm-Stage2.tsv", "\n")
write.table(St2, file =paste(RDATA,"norm-Stage2.tsv", sep="/"), sep="\t", quote=FALSE, row.names=FALSE)
cat("Saving", "norm-Stage3.tsv", "\n")
write.table(St3, file =paste(RDATA,"norm-Stage3.tsv", sep="/"), sep="\t", quote=FALSE, row.names=FALSE)
# cat("Saving", "norm-Stage4.tsv", "\n")
# write.table(St4, file =paste(RDATA,"norm-Stage4.tsv", sep="/"), sep="\t", quote=FALSE, row.names=FALSE)

cat("Saving", "Norm_cpm10_arsyn.RData", "\n")
save(norm.data_cpm10_arsyn, file=paste(RDATA,"Norm_cpm10_arsyn.RData", sep="/"), compress="xz")

cat("End of normalization texting\n")


