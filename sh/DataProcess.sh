#!/bin/bash
echo
echo '##############################################'
echo '########### Preprocessing Data ###############'
echo '##############################################'
echo

#Rscript --vanilla R/01-Rectify_rnaSeq.R Adeno
Rscript --vanilla R/01-Rectify_rnaSeq.R Squamous

echo
echo '###############################################'
echo '######## Pre-Quality Control & Plots ##########'
echo '###############################################'
echo

#Rscript --vanilla R/02-PRE-QC.R Adeno
Rscript --vanilla R/02-PRE-QC.R Squamous

echo
echo '###############################################'
echo '############## Normalization ##################'
echo '###############################################'
echo

#Rscript --vanilla R/03-NORM.R Adeno
Rscript --vanilla R/03-NORM.R Squamous

echo
echo '###############################################'
echo '######## Pre-Quality Control & Plots ##########'
echo '###############################################'
echo

#Rscript --vanilla R/04-POST-QC.R Adeno
#Rscript --vanilla R/04-POST-QC.R Squamous

echo
echo '###############################################'
echo '########### END OF THE PROCESS ################'
echo '###############################################'
echo



