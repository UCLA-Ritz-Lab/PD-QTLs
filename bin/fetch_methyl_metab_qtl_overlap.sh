#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage [peg1cases|peg2cases]"
  exit 1
fi

peg=$1

sql_pd_qtl << END
  select cis.snp_id,cis.allele,count(*) as associated_probes from me_qtls as cis , m_qtls as metab where cis.peg="${peg}" and cis.qtl_type='cis' and
 metab.peg="${peg}" and metab.qtl_type='all' and cis.snp_id=metab.snp_id and cis.allele=metab.allele group by cis.snp_id,cis.allele order by associated_probes desc;
END
