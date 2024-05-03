#!/bin/bash

sql_pd_qtl << END
select trans_peg1.qtl as qtl,trans_peg1.total as trans_peg1_total,trans_peg2.total as trans_peg2_total from qtl_hotspots as trans_peg1 left join qtl_hotspots as trans_peg2 on trans_peg1.qtl_type='trans' and trans_peg1.peg='peg1' and trans_peg2.qtl_type='trans' and  trans_peg2.peg='peg2' and trans_peg1.qtl=trans_peg2.qtl order by trans_peg1_total+trans_peg2_total desc;
#select metab_peg1.qtl as qtl,metab_peg1.total as metab_peg1_total,metab_peg2.total as metab_peg2_total from qtl_hotspots as metab_peg1 left join qtl_hotspots as metab_peg2 on metab_peg1.qtl_type='metab' and metab_peg1.peg='peg1' and metab_peg2.qtl_type='metab' and  metab_peg2.peg='peg2' and metab_peg1.qtl=metab_peg2.qtl order by metab_peg1_total+metab_peg2_total desc;
END
