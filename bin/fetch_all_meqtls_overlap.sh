#!/bin/bash

if [ $# -lt 1 ]
then
  echo "Usage: [1|2]"
  exit 1
fi
peg=$1

sql_pd_qtl << END
select peg1.qtl_type,peg1.snp_id,peg1.allele,peg1.gene,peg1.peg,peg1.pvalue,peg1.FDR,peg1.beta,peg1.beta/peg1.statistic as peg1_se,peg2.peg,peg2.pvalue,peg2.FDR,peg2.beta,peg2.beta/peg2.statistic as peg2_se,meta.zscore as meta_z_score, meta.pvalue as meta_pvalue from me_qtls as peg1,me_qtls as peg2,meta_analysis as meta where peg1.peg='peg1' and peg2.peg='peg2' and peg1.qtl_type=peg2.qtl_type and peg1.snp_id=peg2.snp_id and peg1.allele=peg2.allele and peg1.gene=peg2.gene and peg${peg}.FDR<.05 and meta.qtl_type=peg1.qtl_type and meta.snp_id=peg1.snp_id and meta.gene=peg1.gene order by meta.pvalue;
END
