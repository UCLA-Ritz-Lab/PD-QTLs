#!/usr/bin/env python3
"""
qc_report.py
Summarises per-chromosome VCF QC metrics after pre-processing.
Called by the qc_report rule in the Snakefile.
"""
import os
import sys
import glob
import argparse
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--stats_dir",  required=True)
    p.add_argument("--out_report", required=True)
    return p.parse_args()

def main():
    args = parse_args()

    stat_files = sorted(glob.glob(
        os.path.join(args.stats_dir, "chr*.stats")
    ))

    if not stat_files:
        print(f"No .stats files found in {args.stats_dir}")
        sys.exit(1)

    rows = []
    for f in stat_files:
        chrom = os.path.basename(f).replace(".stats", "").replace("chr", "")
        with open(f) as fh:
            stats = {}
            for line in fh:
                if line.startswith("SN"):
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        stats[parts[1]] = parts[2]
        rows.append({
            "chromosome":  chrom,
            "n_samples":   stats.get("0\tnumber of samples:", "NA"),
            "n_snps":      stats.get("0\tnumber of SNPs:", "NA"),
            "n_indels":    stats.get("0\tnumber of indels:", "NA"),
            "n_multiall":  stats.get("0\tnumber of multiallelic sites:", "NA"),
        })

    df = pd.DataFrame(rows)

    with open(args.out_report, "w") as fh:
        fh.write("=== PEG Imputation Pre-Processing QC Report ===\n\n")
        fh.write("Per-chromosome variant counts after QC:\n\n")
        fh.write(df.to_string(index=False))
        fh.write("\n\n")

        # Totals
        for col in ["n_snps", "n_indels"]:
            try:
                total = sum(int(v.replace(",",""))
                            for v in df[col] if v != "NA")
                fh.write(f"Total {col}: {total:,}\n")
            except Exception:
                pass

        fh.write("\nFiles are ready for upload to:\n")
        fh.write("  Michigan Imputation Server: "
                 "https://imputationserver.sph.umich.edu\n")
        fh.write("  Panel: TOPMed r2 | Build: hg19 | Population: Mixed\n")

    print(f"QC report written to {args.out_report}")
    print(df.to_string(index=False))

if __name__ == "__main__":
    main()
