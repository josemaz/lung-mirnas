library("BiocParallel")
library("NOISeq")

###############################################################################
##Quality Control
###############################################################################
register(MulticoreParam(workers=detectCores()-1, progress=TRUE))

cat("#################\n")
cat("Step 2: QC\n")
cat("#################\n")
DATADIR <- "pipeline/"
RDATA <- paste(DATADIR, "rdata", sep="")
PLOTSDIR <-paste(DATADIR, "plots", sep="")
dir.create(PLOTSDIR)
#PLOTSDIR <-paste(DATADIR, "plots", sep = "")
PREDIR <- paste(PLOTSDIR, "QC_PRE", sep = "/")
dir.create(PREDIR)
w <- 1024 #Resolucion de los plots
h <- 1024  #Resolucion de los plots
p <- 24   #Resolucion de los plots


load(file=paste(RDATA, "RawFull.RData", sep="/"))

##########################################
###Let's keep only the GC & length annotated genes
##########################################
ids<-!is.na(full$Annot$GC) & !is.na(full$Annot$Length)
#ids es un vector de TRUE y FALSE que servirá para quedarse solo con las cuentas de aquellos genes cuya anotacion
#tiene %GC y longitud del gene

full$M<-full$M[ids,] #Aplicamos el filtro creado arriba y quitamos aquellas cuentas que no cumplan con los requisitos
full$Annot<-full$Annot[ids,] #Aplicamos el filtro a las anotaciones y quitamos aquellas anotaciones que no cumplan con los requisitos
cat("Non GC and lenght annotated genes removed.\n")

row.names(full$M)<-full$Annot$EnsemblID #Las cuentas te quedan organizaditas, el gen al que pertenecen, si son de sano o de enfermo y el numero de muestra
row.names(full$Annot)<-full$Annot$EnsemblID #Ahora el numbre de las filas corresponde al ensembleID (lo tenemos 2 veces, una en el nombre de las filas y otra en la columna "esembleID")
row.names(full$Targets)<-full$Targets$ID #El nombre de las filas ahora corresponde al ID (Ctr1, tBa2, etc...)
full$Targets$Group<-factor(substr(full$Targets$ID, start=1, stop=3)) #subtr regresa Ctr y tBa de Ctr20 y tBa21, p.e.
#En targets s crea un acolumna que indica si la muestra es de individuo sano o canceroso

write.table(full$Targets, file=paste(RDATA, "Targets.tsv" , sep="/"), sep="\t", quote=FALSE, row.names=FALSE)
save(full, file=paste(RDATA, "RawFull.RData", sep="/"), compress="xz")
cat("Saving RawFull.RData \n")

##########################################
## EXPLORATORY ANALYSIS (NOISeq package)
##########################################
mydata <- NOISeq::readData(
  data = full$M, 
  length = full$Annot[, c("EnsemblID", "Length")], 
  biotype = full$Annot[, c("EnsemblID", "Type")], 
  chromosome = full$Annot[, c("Chr", "Start", "End")], 
  factors = full$Targets[, "Group",drop=FALSE], 
  gc = full$Annot[, c("EnsemblID", "GC")])

##########################################
## Plots
##########################################
# Biodetection plot. Per group.
mybiodetection <- dat(mydata, type="biodetection", factor="Group", k=0)
png(filename=paste(PREDIR, "01-biodetection.Rd_%02d.png", sep="/"),  width=w, height=h, pointsize=p)
explo.plot(mybiodetection, factor="Group" )
dev.off()
cat("Biodetection plots generated\n")
#Se generan 2 plots (sanos y enfermos) que cotienen el porcentaje de "type" de los archivos 

## Count distribution per biotype. Using count per million, only for one sample
mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
png(filename=paste(PREDIR, "02-countsbio.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")
dev.off()
cat("Counts distribution plot per biotype and one sample generated\n")
#Se genera un plot donde se muestran valores de expresion de los "type"

## Count distribution per sample
png(paste(PREDIR, "02-protein_coding_boxplot.png", sep="/"), width=w*2, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "boxplot")
dev.off()
cat("Counts distribution plot for protein coding and all samples generated\n")
#Se genera un boxplot que muestra valores de expresion de mustras de "protein coding" de sanos y de enfermos

png(paste(PREDIR, "02-protein_coding_barplot.png", sep="/"), width=w*2, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "barplot")
dev.off()
cat("Counts distribution barplot for protein coding biotype and all samples generated\n")
#Se genera un barplot que contiene CPM de genes "protein coding"  de sanos y de enfermos

## Count distribution per Experimental factors
mycountsbio <- dat(mydata, factor = "Group", type = "countsbio")
png(paste(PREDIR, "03-protein_coding_boxplot_group.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "boxplot")
dev.off()
cat("Counts distribution boxplot for protein coding biotype and group generated\n")
#Se genera una grafica con 2 boxplots (sanos y enfermos) que grafican valores de expresion de "protein coding genes"

png(paste(PREDIR, "04-protein_coding_barplot_group.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "barplot")
dev.off()
cat("Counts distribution barplot for protein coding biotype and group generated\n")

##########################################
##Bias
##########################################
## Length bias detection
mylengthbias <- dat(mydata, factor="Group", type="lengthbias")
png(paste(PREDIR, "05-Lengthbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mylengthbias, samples = NULL, toplot = "global")
dev.off()
cat("Lenght bias plot generated\n")

##GC bias
mygcbias <- dat(mydata, factor = "Group", type="GCbias")
png(paste(PREDIR, "06-GCbias.png", sep="/"), width=w, height=h, pointsize=p)
explo.plot(mygcbias, samples = NULL, toplot = "global")
dev.off()
cat("GC bias plot generated\n")

## RNA composition
#mycomp <- dat(mydata, type="cd")
# png(paste(PREDIR, "07-RNAComposition.png", sep="/"), width=w, height=h, pointsize=p)
# explo.plot(mycomp, samples=1:12)
# dev.off()
# cat("RNA composition plot generated\n")

#########################
# Quality Control Report
#########################
# A complete pdf report can be obtained using this function.
#¿No funciona con mas de 2 grupos?
QCreport(mydata, factor="Group", file=paste(PLOTSDIR, "QCReport.pdf", sep="/"))


#############################
## PCA Analysis with NOISeq
##########################################
pca.dat <- dat(mydata, type = "PCA", logtransf = F)
pca.results <- pca.dat@dat$result

## Variance explained by each component
pdf(file=paste(PREDIR, "08-PCAVariance_raw.pdf", sep="/"), width = 4*2, height = 4*2)
barplot(pca.results$var.exp[,1], xlab = "PC", ylab = "Explained variance")
dev.off()
cat("PCA variance raw plot generated.\n")

## Loading plot
pdf(file=paste(PREDIR, "09-PCALoading_raw.pdf", sep="/"), width = 4*2, height = 4*2)
plot(pca.results$loadings[,1:2], col = 1, pch = 20, cex = 0.5,
     xlab = paste("PC 1 ", round(pca.results$var.exp[1,1]*100,0), "%", sep = ""),
     ylab = paste("PC 2 ", round(pca.results$var.exp[2,1]*100,0), "%", sep = ""),
     main = "PCA loadings",
     xlim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1),
     ylim = range(pca.results$loadings[,1:2]) + 0.02*diff(range(pca.results$loadings[,1:2]))*c(-1,1))  
dev.off()
cat("PCA loading raw plot generated.\n")

## Score plot
mycol <- as.character(full$Targets$Group)
mycol[mycol == 'Ctr'] <- "red2"
mycol[mycol == 'tBa'] <- "blue1"
# mycol[mycol == 's4-'] <- "cyan2"

pdf(file=paste(PREDIR, "10-PCAScore_raw.pdf", sep="/"), width = 5*2, height = 5)
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
legend("topright", c("Ctr","tBa"), 
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
legend("topright", c("Ctr", "tBA"), 
       # col = c("black", "red2", "blue1", "green2", "cyan2"), ncol = 2, pch = 1)
       col = c("red2", "blue1"), ncol = 2, pch = 1)
dev.off()
cat("PCA scores raw plot generated.\n")

##########################################################################
## GREAT JOB!!! YOU MADE IT TILL THE END!!!!
###########################################################################
cat("#################\n")
cat("End of Step 2: QC\n")
cat("#################\n")

