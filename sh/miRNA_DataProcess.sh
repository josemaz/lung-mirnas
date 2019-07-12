#!/bin/bash
echo
echo '####################################################'
echo '########### Preprocessing miRNA Data ###############'
echo '####################################################'
echo

echo '--------------------- Adeno ----------------------'
Rscript --vanilla R/miRNA-Normalize.R Adeno
echo
echo '------------------- Squamous ---------------------'
Rscript --vanilla R/miRNA-Normalize.R Squamous
echo

echo
echo '#####################################################'
echo '############ END OF miRNA PROCESSING ################'
echo '#####################################################'
echo
