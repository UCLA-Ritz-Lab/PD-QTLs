#!/bin/bash

# this code will fetch all meQTLs that are found in both PEG1 and PEG2 cases but not in controls

sql_pd_qtl << END
select peg1cases.qtl_type,concat(peg1cases.snp_id,',',peg1cases.gene) as 'snp,gene',snpinfo.chrom,snpinfo.pos,snpinfo.ref_allele,snpinfo.alt_allele,peg1cases.peg,peg1cases.pvalue,peg1cases.FDR,peg1cases.beta,peg1cases.beta/peg1cases.statistic as se,peg2cases.peg,peg2cases.pvalue,peg2cases.FDR,peg2cases.beta,peg2cases.beta/peg2cases.statistic as se from me_qtls as peg1cases inner join me_qtls as peg2cases on (peg2cases.peg='peg2cases' and peg1cases.snp_id=peg2cases.snp_id and peg1cases.gene=peg2cases.gene) inner join snpinfo on (snpinfo.snpid2=peg1cases.snp_id) left join me_qtls as peg1controls on (peg1controls.peg='peg1controls' and peg1cases.snp_id is null) where peg1cases.peg='peg1cases' and peg1cases.FDR<.05 and peg2cases.FDR<.05;
END
