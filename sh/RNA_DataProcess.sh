#!/bin/bash
echo
echo '##################################################'
echo '########### Preprocessing RNA Data ###############'
echo '##################################################'
echo

echo '--------------------- Adeno ----------------------'
Rscript --vanilla R/01-Rectify_rnaSeq.R Adeno
echo
echo '------------------- Squamous ---------------------'
Rscript --vanilla R/01-Rectify_rnaSeq.R Squamous
echo

echo
echo '###############################################'
echo '######## Pre-Quality Control & Plots ##########'
echo '###############################################'
echo

echo '--------------------- Adeno ----------------------'
Rscript --vanilla R/02-PRE-QC.R Adeno
echo
echo '------------------- Squamous ---------------------'
Rscript --vanilla R/02-PRE-QC.R Squamous
echo

echo
echo '###############################################'
echo '############## Normalization ##################'
echo '###############################################'
echo

echo '--------------------- Adeno ----------------------'
Rscript --vanilla R/03-NORM.R Adeno
echo
echo '------------------- Squamous ---------------------'
Rscript --vanilla R/03-NORM.R Squamous
echo

echo
echo '###############################################'
echo '######## Pre-Quality Control & Plots ##########'
echo '###############################################'
echo

echo '--------------------- Adeno ----------------------'
Rscript --vanilla R/04-POST-QC.R Adeno
echo
echo '------------------- Squamous ---------------------'
Rscript --vanilla R/04-POST-QC.R Squamous
echo

echo
echo '###############################################'
echo '########### END OF RNA PROCESS ################'
echo '###############################################'
echo



