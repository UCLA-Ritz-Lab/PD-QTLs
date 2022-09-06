#!/bin/bash

sql_pd_qtl << END
select chrom as CHR,pos as BP,snp_id as SNP,p_value as P from peg1_gwas as a order by p_value;
END

