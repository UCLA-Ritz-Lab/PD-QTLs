#!/bin/bash

sql_pd_qtl << END
select IlmnID,me_qtls.snp_id as snp_id,peg1_gwas.pos as position,illumina_annotation.chrom,mapinfo-1e6,mapinfo-1,me_qtls.beta as me_qtl_beta,power(me_qtls.beta/me_qtls.statistic,2) as me_qtl_var, log(peg1_gwas.odd_ratio) as gwas_beta,power(beta_se,2) as gwas_var from illumina_annotation,me_qtls,peg1_gwas where me_qtls.gene=IlmnID and peg1_gwas.snp_id=me_qtls.snp_id and peg='peg1' and me_qtl_type='cis' and test='ADD' order by IlmnID;
END
