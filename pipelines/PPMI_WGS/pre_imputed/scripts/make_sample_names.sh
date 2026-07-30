#!/bin/bash

#Prepare sample manifest

cd ../../../Downloads/PPMI_WGS/pre_imputed
ls PPMISI*.realigned*gz|cut -f1 -d'.'|sed 's/$/.variant2/' > ../../../results/PPMI_WGS/pre_imputed/samples_snps.txt
ls PPMISI*.realigned*gz|cut -f1 -d'.'|sed 's/$/.variant/' > ../../../results/PPMI_WGS/pre_imputed/samples_indels.txt
cd ../../../pipelines/PPMI_WGS/pre_imputed
