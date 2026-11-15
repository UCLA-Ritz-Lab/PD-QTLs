# best_allele.awk
# Reads a PLINK2 .afreq and selects, for each chr:pos position, the single
# variant record whose ALT has the highest cohort-wide frequency.
#
# Input : #CHROM  ID(chr:pos:ref:alt)  REF  ALT  ALT_FREQS  OBS_CT
# Output: <full variant ID>\t<ALT allele of that same record>
#
# The output doubles as the --export-allele file and (via its first column)
# the --extract list, so only the selected record survives to export.
#
# Why one record per position rather than mapping every record at a position
# to the position's best allele: after multiallelic splitting, each record is
# biallelic, so a record with alleles C/A cannot count allele G. PLINK2 accepts
# such a request without warning and emits all-zero dosage, and the downstream
# chr:pos dedup then keeps whichever record came first — sometimes the zeroed
# one. Selecting the record here means the exported allele is always genuinely
# present, and no two surviving records share a chr:pos.
#
# Records at the same position are adjacent (the .afreq follows pvar/genomic
# order), so grouping is done in a single pass with O(1) memory rather than
# hashing all ~21M variant IDs.
/^#/ { next }
{
    split($2, p, ":")
    key = p[1] ":" p[2]

    if (key != cur) {
        if (cur != "") print best_id "\t" best_alt
        cur      = key
        best_id  = $2
        best_alt = $4
        best_frq = $5 + 0
        next
    }
    # Strictly greater, so the first record wins ties — deterministic given
    # the fixed input order.
    if ($5 + 0 > best_frq) {
        best_id  = $2
        best_alt = $4
        best_frq = $5 + 0
    }
}
END {
    if (cur != "") print best_id "\t" best_alt
}
