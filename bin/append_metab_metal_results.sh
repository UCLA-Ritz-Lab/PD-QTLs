#!/bin/bash

sql_pd_qtl  << END
END
echo "load data infile '/var/analysis/PD-QTLs/results/met_qtls/metal.import' into table meta_analysis;" | sql_pd_qtl
