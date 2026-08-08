#!/usr/bin/env bash
# Build a cohort's sample keep-list for the ancestry PCA.
#
# The keep-list has to be in the VCF's *own* sample-ID namespace, because it is
# fed straight to `bcftools view -S` by rule cohort_pgen. That is trivial when a
# cohort's covariate file is keyed by the same IDs as its genotypes (PPMI), but
# PEG keys the two differently: the VCF carries GWAS_IDs (CRG_*/CRG2_*) while
# every PEG covariate file is keyed by Pegid. Intersecting those two namespaces
# directly yields the empty set. Such cohorts declare a linkage file and we join
# through it.
#
# Two outputs, deliberately:
#   --keep   VCF-space IDs, one per line, for `bcftools view -S`
#   --id-map vcf_id,cohort_id CSV, so downstream stages that read covariates
#            (which are cohort-space) can line them up against PCA output
#            (which is VCF-space). Without this the self-report join silently
#            matches nothing and every sample lands in the "Unknown" bucket.
#
# CSV parsing is naive on purpose (plain `-F,`): none of these inputs quote
# embedded commas in their ID columns. If that ever changes this needs a real
# CSV reader.
set -euo pipefail

cohort= vcf= cov= cov_col=
linkage= link_from= link_to= link_dup_col=
keep= idmap= report=

die() { echo "ERROR: $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cohort)       cohort=$2;       shift 2 ;;
    --vcf)          vcf=$2;          shift 2 ;;
    --cov)          cov=$2;          shift 2 ;;
    --cov-column)   cov_col=$2;      shift 2 ;;
    --linkage)      linkage=$2;      shift 2 ;;
    --link-from)    link_from=$2;    shift 2 ;;
    --link-to)      link_to=$2;      shift 2 ;;
    --link-dup-column) link_dup_col=$2; shift 2 ;;
    --keep)         keep=$2;         shift 2 ;;
    --id-map)       idmap=$2;        shift 2 ;;
    --report)       report=$2;       shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done

for req in cohort vcf cov cov_col keep idmap report; do
  [[ -n "${!req}" ]] || die "missing required --${req//_/-}"
done

mkdir -p "$(dirname "$keep")" "$(dirname "$idmap")" "$(dirname "$report")"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# ── Helpers ───────────────────────────────────────────────────────────────────

# csv_col <file> <col-name>            -> one column, header dropped
# csv_cols <file> <col-a> <col-b> [<dup-col>] -> "a<TAB>b", dropping dup-col != 0
csv_col() {
  awk -F, -v want="$2" -v f="$1" '
    NR == 1 {
      for (i = 1; i <= NF; i++) { h = $i; gsub(/^"|"$/, "", h); if (h == want) c = i }
      if (!c) { print "ERROR: no column \"" want "\" in " f > "/dev/stderr"; exit 1 }
      next
    }
    { v = $c; gsub(/^"|"$/, "", v); if (v != "") print v }
  ' "$1"
}

csv_cols() {
  awk -F, -v wa="$2" -v wb="$3" -v wd="${4:-}" -v f="$1" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        h = $i; gsub(/^"|"$/, "", h)
        if (h == wa) ca = i; if (h == wb) cb = i; if (wd != "" && h == wd) cd = i
      }
      if (!ca) { print "ERROR: no column \"" wa "\" in " f > "/dev/stderr"; exit 1 }
      if (!cb) { print "ERROR: no column \"" wb "\" in " f > "/dev/stderr"; exit 1 }
      if (wd != "" && !cd) { print "ERROR: no column \"" wd "\" in " f > "/dev/stderr"; exit 1 }
      next
    }
    {
      a = $ca; gsub(/^"|"$/, "", a)
      b = $cb; gsub(/^"|"$/, "", b)
      if (a == "" || b == "") next
      if (cd) { d = $cd; gsub(/^"|"$/, "", d); if (d != "0") next }
      print a "\t" b
    }
  ' "$1"
}

# ── Gather IDs ────────────────────────────────────────────────────────────────

bcftools query -l "$vcf" | sort -u > "$tmp/vcf_ids"
csv_col "$cov" "$cov_col" | sort -u > "$tmp/cov_ids"

if [[ -n "$linkage" ]]; then
  [[ -n "$link_from" && -n "$link_to" ]] \
    || die "--linkage requires --link-from and --link-to"

  # cohort-space -> VCF-space, restricted to the cohort's covariate IDs
  csv_cols "$linkage" "$link_from" "$link_to" "$link_dup_col" | sort -u > "$tmp/link"
  join -t$'\t' "$tmp/cov_ids" "$tmp/link" | sort -u > "$tmp/cov_pairs"

  # Keep only pairs whose VCF-space ID is genuinely in the VCF.
  awk -F'\t' 'NR == FNR { v[$1]; next } ($2 in v)' \
    "$tmp/vcf_ids" "$tmp/cov_pairs" | sort -u > "$tmp/pairs"
else
  # Covariates are already in VCF space: the pair is the ID with itself.
  comm -12 "$tmp/vcf_ids" "$tmp/cov_ids" | awk '{ print $1 "\t" $1 }' > "$tmp/pairs"
fi

cut -f2 "$tmp/pairs" | sort -u > "$keep"

{
  echo "vcf_id,cohort_id"
  # vcf_id first: that is the key downstream stages join PCA output on.
  awk -F'\t' '{ print $2 "," $1 }' "$tmp/pairs" | sort -u
} > "$idmap"

# ── Report ────────────────────────────────────────────────────────────────────

n_vcf=$(wc -l < "$tmp/vcf_ids")
n_cov=$(wc -l < "$tmp/cov_ids")
n_keep=$(wc -l < "$keep")

{
  echo "Sample selection: $cohort"
  echo "  samples in VCF          : $n_vcf"
  echo "  samples in covariates   : $n_cov  (column: $cov_col)"
  if [[ -n "$linkage" ]]; then
    echo "  linkage                 : $linkage"
    echo "    $link_from -> $link_to$([[ -n $link_dup_col ]] && echo ", dropping $link_dup_col != 0")"
    echo "    usable linkage rows   : $(wc -l < "$tmp/link")"
  else
    echo "  linkage                 : none (covariate IDs are VCF sample IDs)"
  fi
  echo "  kept (VCF-space IDs)    : $n_keep"
  echo
  echo "  covariate IDs with no VCF genotypes:"
  comm -23 "$tmp/cov_ids" <(cut -f1 "$tmp/pairs" | sort -u) | sed 's/^/    /'
} > "$report"

[[ "$n_keep" -gt 0 ]] || die "no samples survived the intersection — see $report"

# A cohort-space ID mapping to several VCF samples would enter the PCA more than
# once and drag the affected region of the plot toward itself. Fail loudly.
dupes=$(cut -f1 "$tmp/pairs" | sort | uniq -d)
if [[ -n "$dupes" ]]; then
  echo "ERROR: $cohort: these $cov_col values map to >1 VCF sample:" >&2
  echo "$dupes" | sed 's/^/  /' >&2
  die "refusing to build a keep-list with duplicated individuals"
fi
