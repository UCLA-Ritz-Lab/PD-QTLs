# psam_to_keep.awk
# Emits a PLINK2 --keep file (FID IID) covering every sample in a .psam.
#
# Handles both psam layouts seen in this project:
#   Two-col (#FID + IID):  imputed pfiles — FID=IID=PPMISI{PATNO}
#   Single-col (#IID):     old WGS pfiles — FID assigned as 0 by PLINK2
#
# Note there is no sample-level filtering here: the imputed psam contains only
# real genotyped samples. The pre-imputation WGS psam additionally held
# non-genotype ".variant" rows, which is why this step used to grep for the
# ".variant2" suffix — the imputation server dropped that suffix, so the grep
# matched nothing and silently produced an empty keep file.
NR==1 && /^#/ {
    for (i = 1; i <= NF; i++) {
        h = $i
        sub(/^#/, "", h)
        if (h == "FID") fid = i
        if (h == "IID") iid = i
    }
    if (!iid) iid = 1
    next
}
{
    if (fid) print $fid, $iid
    else     print 0, $iid
}
