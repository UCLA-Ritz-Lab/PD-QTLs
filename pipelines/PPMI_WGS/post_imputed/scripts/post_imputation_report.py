#!/usr/bin/env python3
"""
post_imputation_report.py

Summarises HRC imputation quality using per-chromosome .info.gz files
returned by the Michigan Imputation Server, plus final variant count
from the merged filtered VCF.

Usage:
    python post_imputation_report.py \
        --info_dir results/unzipped \
        --final_vcf results/peg_imputed_filtered.dose.vcf.gz \
        --rsq_min 0.3 \
        --out results/qc/post_imputation_report.txt
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--info_dir",  required=True, help="Directory containing chr{N}.info.gz")
    p.add_argument("--final_vcf", required=True, help="Merged filtered VCF")
    p.add_argument("--mis_qc",    default="../../../results/PPMI_WGS/imputation_results/qc_report.txt", help="MIS qc_report.txt path")
    p.add_argument("--rsq_min",     type=float, default=0.3)
    p.add_argument("--out",       required=True, help="Output report path")
    return p.parse_args()



def parse_mis_qc_report(path):
    """Parse the MIS pre-imputation qc_report.txt into a dict of key/value stats."""
    stats = {}
    if not os.path.exists(path):
        return stats
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            # Strip HTML bold tags
            line = line.replace("<b>", "").replace("</b>", "")
            # Skip directive lines and section markers
            if line.startswith("::") or not line:
                continue
            # Lines like "Match: 422,482" or "Reference Overlap: 100.00 %"
            if ":" in line:
                key, _, val = line.partition(":")
                stats[key.strip()] = val.strip()
    return stats

def parse_info_file(path):
    """Parse a minimac4 .info.gz file into a DataFrame.
    Despite the .info.gz extension, MIS returns a BGZF-compressed VCF file.
    Imputation quality metrics (R2, MAF, AvgCall) are in the INFO column,
    e.g. AF=0.123;MAF=0.123;R2=0.987;AvgCall=0.998;Imputed
    """
    import gzip
    import re

    records = []
    info_pat = re.compile(r'(?:^|;)(\w+)=([^;]+)')

    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, snp_id, ref, alt, qual, filt, info_str = parts[:8]
            info_dict = dict(info_pat.findall(info_str))
            records.append({
                "SNP":    snp_id,
                "CHROM":  chrom,
                "POS":    pos,
                "REF":    ref,
                "ALT":    alt,
                "Rsq":    info_dict.get("R2"),
                "MAF":    info_dict.get("MAF"),
                "AvgCS":  info_dict.get("AVG_CS"),
                "Imputed": "IMPUTED" in info_str,
            })

    df = pd.DataFrame(records)
    for col in ["Rsq", "MAF", "AvgCS"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def count_vcf_variants(vcf_path):
    """Count SNPs in a VCF using bcftools stats."""
    result = subprocess.run(
        ["bcftools", "stats", vcf_path],
        capture_output=True, text=True, check=True
    )
    for line in result.stdout.splitlines():
        if line.startswith("SN") and "number of SNPs" in line:
            return int(line.split("\t")[3])
    return None


def main():
    args = parse_args()
    info_dir = Path(args.info_dir)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    chroms = list(range(1, 23))
    per_chrom = []

    print("Parsing .info.gz files...", file=sys.stderr)
    for chrom in chroms:
        info_file = info_dir / f"chr{chrom}.info.gz"
        if not info_file.exists():
            print(f"  WARNING: {info_file} not found, skipping", file=sys.stderr)
            continue

        df = parse_info_file(info_file)
        total      = len(df)
        pass_rsq   = (df["Rsq"] >= args.rsq_min).sum()
        well_imp   = (df["Rsq"] >= 0.8).sum()
        median_rsq = df["Rsq"].median()
        mean_rsq   = df["Rsq"].mean()

        per_chrom.append({
            "chrom":         chrom,
            "total_variants": total,
            "pass_rsq":      pass_rsq,
            "well_imputed":  well_imp,
            "median_rsq":    round(median_rsq, 4),
            "mean_rsq":      round(mean_rsq, 4),
            "pct_pass":      round(100 * pass_rsq / total, 1) if total > 0 else 0
        })

    summary_df = pd.DataFrame(per_chrom)

    # Overall stats
    total_all    = summary_df["total_variants"].sum()
    pass_all     = summary_df["pass_rsq"].sum()
    well_all     = summary_df["well_imputed"].sum()

    # Count final merged VCF
    print("Counting final VCF variants...", file=sys.stderr)
    final_count = count_vcf_variants(args.final_vcf)

    lines = [
        "=" * 60,
        "PPMI POST-IMPUTATION QC REPORT",
        "Reference Panel : HRC r1.1 2016 (GRCh37/hg19)",
        f"Rsq filter      : >= {args.rsq_min}",
        "=" * 60,
        "",
        "── MIS pre-imputation QC ───────────────────────────────────",
    ]

    mis_qc = parse_mis_qc_report(args.mis_qc)
    if mis_qc:
        for k, v in mis_qc.items():
            lines.append(f"  {k:<40}: {v}")
    else:
        lines.append("  (qc_report.txt not found)")

    lines += [
        "",
        "── Per-chromosome summary ──────────────────────────────────",
        f"{'Chr':<6} {'Total':>10} {'Pass Rsq':>10} {'Well-imp':>10} {'%Pass':>7} {'Median Rsq':>12}",
        "-" * 60,
    ]

    for row in per_chrom:
        lines.append(
            f"{row['chrom']:<6} {row['total_variants']:>10,} "
            f"{row['pass_rsq']:>10,} {row['well_imputed']:>10,} "
            f"{row['pct_pass']:>6.1f}% {row['median_rsq']:>12.4f}"
        )

    lines += [
        "-" * 60,
        f"{'TOTAL':<6} {total_all:>10,} {pass_all:>10,} {well_all:>10,} "
        f"{100*pass_all/total_all:>6.1f}%",
        "",
        "── Overall statistics ──────────────────────────────────────",
        f"Total variants (pre-filter)   : {total_all:,}",
        f"Variants passing Rsq >= {args.rsq_min}    : {pass_all:,}  ({100*pass_all/total_all:.1f}%)",
        f"Well-imputed variants (Rsq>=0.8): {well_all:,}  ({100*well_all/total_all:.1f}%)",
        f"Final merged VCF variant count: {final_count:,}" if final_count else "",
        "",
        "── Notes ───────────────────────────────────────────────────",
        "Well-imputed (Rsq >= 0.8) variants are suitable for",
        "primary mQTL analysis. Variants with 0.3 <= Rsq < 0.8",
        "may be retained but should be interpreted with caution.",
        "=" * 60,
    ]

    report = "\n".join(lines)
    print(report)

    with open(out_path, "w") as fh:
        fh.write(report + "\n")

    print(f"\nReport written to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
