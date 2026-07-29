BEGIN { OFS="\t" }
NR==1 {
    printf "snp_id"
    for (i=7; i<=NF; i++) {
        id = $i
        sub(/^[^_]+_/, "", id)
        printf "\t%s", id
    }
    printf "\n"
}
NR>1 {
    id = $2
    sub(/_[ACGT]+$/, "", id)
    printf "%s", id
    for (i=7; i<=NF; i++) printf "\t%s", $i
    printf "\n"
}
