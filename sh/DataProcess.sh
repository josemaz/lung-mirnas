#!/bin/bash
echo
echo '##############################################'
echo '########### Preprocessing Data ###############'
echo '##############################################'
echo

Rscript --vanilla R/01-Rectify_rnaSeq.R Adeno
Rscript --vanilla R/01-Rectify_rnaSeq.R Scamous

echo
echo '###############################################'
echo '######## Pre-Quality Control & Plots ##########'
echo '###############################################'
echo
Rscript --vanilla R/02-PRE-QC.R Adeno
Rscript --vanilla R/02-PRE-QC.R Scamous
