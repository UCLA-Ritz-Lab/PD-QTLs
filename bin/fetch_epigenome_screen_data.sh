#!/bin/bash

echo 'select b.peg_id,PDstudyParkinsonsDisease from datMethPEG1t_sample as a,peg1_covariates as b where b.sample_id=a.sample_id order by a.seq'| sql_pd_qtl > peg1subjects_disease.txt
