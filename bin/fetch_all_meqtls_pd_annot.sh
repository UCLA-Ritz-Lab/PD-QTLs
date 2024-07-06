#!/bin/bash

#if [ $# -lt 1 ]
#then
#  echo "Usage: [1|2]"
#  exit 1
#fi
#peg=$1

sql_pd_qtl << END
#select peg1.me_qtl_type,peg1.snp_id,peg1.allele,peg1.gene,peg1.peg,peg1.pvalue,peg1.FDR,peg1.beta,power(peg1.beta/peg1.statistic,2) as betavar, pdsnps.beta as pd_beta,power(pdsnps.se,2) as pd_betavar,pdsnps.pvalue,pdsnps.gene from me_qtls as peg1 inner join pd_snps_nall as  pdsnps on (peg1.snp_id=pdsnps.snpid ) where peg1.peg='peg1' and peg${peg}.FDR<.05 order by peg${peg}.FDR;
select meta.snp_id,meta.gene as probe,meta.zscore,meta.pvalue,pdsnps.beta as pd_beta,pdsnps.se as pd_beta_se,pdsnps.pvalue,pdsnps.gene,pdsnps_suppl.beta_all_studies as pd_beta,pdsnps_suppl.se_all_studies as pd_beta_se,pdsnps_suppl.pvalue_all_studies from snpinfo inner join meta_analysis as meta on (snpinfo.snpid2=meta.snp_id) inner join pd_snps_nall as pdsnps on (snpinfo.snpid=pdsnps.snpid) left join pd_snps_nall_suppl as pdsnps_suppl on (snpinfo.snpid=pdsnps_suppl.snpid) where qtl_type='cis' order by meta.pvalue;
END
