select d.gwas_id,d.gwas_id from datMethPEG1t_sample as a,peg1_covariates
as b,peg_id_gwas_id_mapping as c,gwas_subjects as d where b.sample_id=a.sample_id and c.peg_id=b.peg_id and c.dup<>1 and d.gwas_id=c.gwas_id order by d.gwas_id;
