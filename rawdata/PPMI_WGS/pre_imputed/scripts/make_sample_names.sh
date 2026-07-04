#!/bin/bash

#Prepare sample manifest

cd Downloads
ls PPMISI*.realigned*gz|cut -f1 -d'.'|sed 's/$/.variant2/' > ../samples_snps.txt
ls PPMISI*.realigned*gz|cut -f1 -d'.'|sed 's/$/.variant/' > ../samples_indels.txt
cd ..
