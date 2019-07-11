###############################################################################
##Get the Tissue Type (Adenocarcinoma or Squamous Cells carcinoma)  
###############################################################################
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Execution: Rscript --vanilla R/04-POST-QC.R (tissue type); tissue type: Adeno or Squamous", call.=FALSE)
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
library("ggplot2")
library("reshape2")
###############################################################################
RDATA <- paste("Data",Tissue, "rdata", "RNA", sep = "/")
PLOTSDIR <-paste("Plots", Tissue, "RNA", sep="/")
POSTDIR <- paste(PLOTSDIR, "QC_POST", sep = "/")
dir.create(POSTDIR)
w <- 1024
h <- 1024
p <- 24


load(file=paste(RDATA, "Norm_cpm10_arsyn.RData", sep="/"))

##########################################
## EXPLORATORY ANALYSIS (NOISeq package)
##########################################
mydata <- NOISeq::readData(
  data = norm.data_cpm10_arsyn$M, 
  length = norm.data_cpm10_arsyn$Annot[, c("EnsemblID", "Length")], 
  biotype = norm.data_cpm10_arsyn$Annot[, c("EnsemblID", "Type")], 
  chromosome = norm.data_cpm10_arsyn$Annot[, c("Chr", "Start", "End")], 
  factors = norm.data_cpm10_arsyn$Targets[, "Group",drop=FALSE], 
  gc = norm.data_cpm10_arsyn$Annot[, c("EnsemblID", "GC")])

### Length bias
mylengthbias <- dat(mydata, factor="Group", norm = TRUE, type="lengthbias")
png(paste(POSTDIR,"01-Lenghtbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mylengthbias, samples = NULL, toplot = "global")
dev.off()
cat("Lenght bias plot generated\n")

## GC Bias
mygcbias <- dat(mydata, factor = "Group", norm = TRUE, type ="GCbias")
png(paste(POSTDIR,"02-GCbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mygcbias, samples = NULL, toplot = "global")
dev.off()
cat("GC bias plot generated\n")

#RNA Composition
# myrnacomp <- dat(mydata, norm = TRUE, type="cd")
# png(paste(POSTDIR,"03-RNAComposition.png", sep="/"), width=w, height=h, pointsize=p)
# explo.plot(myrnacomp, samples = 1:12)
# dev.off()
# cat("RNA composition plot generated\n")

mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
png(paste(POSTDIR,"04-protein_coding_boxplot.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")
dev.off()
cat("Counts distribution plot for protein coding and all samples generated\n")


#############################
## PCA Analysis with NOISeq
##########################################

pca.dat <- dat(mydata, type = "PCA", logtransf = F)
pca.results <- pca.dat@dat$result

## Variance explained by each component
pdf(file=paste(POSTDIR, "05-PCAVariance_raw.pdf", sep="/"), width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance")
dev.off()
cat("PCA variance raw plot generated.\n")

## Loading plot
pdf(file=paste(POSTDIR, "06-PCALoading_raw.pdf", sep="/"), width = 4*2, height = 4*2)
plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA loadings",
     xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
     ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
dev.off()
cat("PCA loading raw plot generated.\n")

## Score plot
mycol <- as.character(norm.data_cpm10_arsyn$Targets$Group)
mycol[mycol == normalTissue] <- "red2"
mycol[mycol == cancerTissue] <- "blue1"

# mycol[mycol == 's4-'] <- "cyan2"

pdf(file=paste(POSTDIR, "07-PCAScore_raw.pdf", sep="/"), width = 5*2, height = 5)
par(mfrow = c(1,2))

# PC1 & PC2
rango <- diff(range(pca.results$scores[,1:2]))
plot(pca.results$scores[,1:2], col = "white",
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA scores",
     xlim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1),
     ylim = range(pca.results$scores[,1:2]) + 0.02*rango*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,2], col = mycol, cex = 1.5)  
# legend("topright", c("Ctr", "s1-", "s2-", "s3-", "s4-"), 
legend("topright", c(normalTissue, cancerTissue), 
       # col = c("black", "red2", "blue1", "green2", "cyan2"), ncol = 2, pch = 1)
       col = c("red2", "blue1"), ncol = 2, pch = 1)

# PC1 & PC3
rango2 = diff(range(pca.results$scores[,c(1,3)]))
plot(pca.results$scores[,c(1,3)], col = "white",
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 3 ", round(pca.results$var.exp[3,1]*100,0), "%", sep = ""),
     main = "PCA scores",
     xlim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1),
     ylim = range(pca.results$scores[,c(1,3)]) + 0.02*rango2*c(-1,1))
points(pca.results$scores[,1], pca.results$scores[,3], col = mycol, cex = 1.5)
# legend("topright", c("Ctr", "s1-", "s2-", "s3-", "s4-"), 
legend("topright", c(normalTissue, cancerTissue), 
       # col = c("black", "red2", "blue1", "green2", "cyan2"), ncol = 2, pch = 1)
       col = c("red2", "blue1"), ncol = 2, pch = 1)
dev.off()
cat("PCA scores raw plot generated.\n")

#############################
##########################################################################
## GREAT JOB!!! YOU MADE IT TILL THE END!!!!
###########################################################################
cat("####################################\n")
cat("End of Step 4: Post Quality Control\n")
cat("####################################\n")
