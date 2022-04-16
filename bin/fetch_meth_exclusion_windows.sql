drop table if exists meth_exclusion;
create temporary table meth_exclusion(IlmnID varchar(255),start bigint,end bigint,primary key(IlmnID));
insert into meth_exclusion select IlmnID,mapinfo,mapinfo from illumina_annotation where chrom is not null and chrom<>'X' and chrom<>'Y';
update meth_exclusion set start=(start-20000);
update meth_exclusion set end=(end+20000);
update meth_exclusion set start=1 where start<1;
select chrom,start,end,a.IlmnID from illumina_annotation as a,meth_exclusion as b where a.IlmnID=b.IlmnID order by chrom,start,end,a.IlmnID;


