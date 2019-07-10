#!/bin/bash
echo
echo '####################################################'
echo '########### Preprocessing miRNA Data ###############'
echo '####################################################'
echo

Rscript --vanilla R/miRNA-Normaliza.R Adeno
Rscript --vanilla R/miRNA-Normalize.R Squamous

echo
echo '#####################################################'
echo '############ END OF miRNA PROCESSING ################'
echo '#####################################################'
echo
