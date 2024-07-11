#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage [peg1cases|peg2cases]"
  exit 1
fi

peg=$1

sql_pd_qtl << END
  select cis.snp_id as common_eqtl,cis.gene as cis_cpg,cis.statistic as mqtl_statistic,cis.pvalue as mqtl_pvalue,cis.FDR as mqtl_FDR,cis.beta as mqtl_beta,metab.gene as metabolite,metab.statistic as met_statistic,metab.pvalue as met_pvalue,metab.FDR as met_fdr,metab.beta as met_beta from me_qtls as cis, met_qtls as metab where cis.peg="${peg}" and cis.qtl_type='cis' and metab.peg="${peg}" and metab.qtl_type='all' and cis.snp_id=metab.snp_id order by cis.FDR desc;
END
