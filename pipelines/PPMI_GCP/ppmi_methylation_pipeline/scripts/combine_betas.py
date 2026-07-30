#!/usr/bin/env python3
"""
combine_betas.py
Combines per-batch beta matrices (from R/SeSAMe) into a single
probe x sample matrix. Handles:
  - probe set intersection across batches (retains only probes
    passing QC in all batches)
  - duplicate sample detection (same SubjectID_T{visit} in two batches)
  - outputs final matrix as both .csv.gz and .parquet for downstream use
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd

from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--betas_dir",     required=True,
                   help="Directory containing per-batch .betas.rds files "
                        "(read as .csv.gz — see Snakefile conversion step)")
    p.add_argument("--qc_dir",        required=True,
                   help="Directory containing per-batch .qc.csv files")
    p.add_argument("--out_matrix",    required=True,
                   help="Output path for combined beta matrix (.csv.gz)")
    p.add_argument("--out_qc",        required=True,
                   help="Output path for combined QC summary (.csv)")
    p.add_argument("--out_probes",    required=True,
                   help="Output path for final retained probe list (.txt)")
    p.add_argument("--min_probe_detection", type=float, default=0.90,
                   help="Min fraction of samples a probe must pass in to retain")
    p.add_argument("--log",           default="combine_betas.log")
    return p.parse_args()


def log(msg, logfile):
    import datetime
    line = f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    print(line, file=logfile, flush=True)


def main():
    args = parse_args()

    with open(args.log, "w") as logfile:

        # ── Discover batch beta CSV files ─────────────────────────────────────
        # Snakefile converts .rds -> .csv.gz before calling this script
        beta_files = sorted(glob.glob(
            os.path.join(args.betas_dir, "*.betas.csv.gz")
        ))
        if not beta_files:
            log("ERROR: No .betas.csv.gz files found in " + args.betas_dir,
                logfile)
            sys.exit(1)

        log(f"Found {len(beta_files)} batch beta files", logfile)

        # ── Load all batch QC files ───────────────────────────────────────────
        qc_files = sorted(glob.glob(os.path.join(args.qc_dir, "*.qc.csv")))
        qc_all = pd.concat(
            [pd.read_csv(f) for f in qc_files],
            ignore_index=True
        )
        log(f"Total samples across all batches: {len(qc_all)}", logfile)

        # Check for duplicate sample column names across batches
        dup_samples = qc_all[
            qc_all.duplicated(subset="sample_col_name", keep=False)
        ]
        if len(dup_samples) > 0:
            log(f"WARNING: {len(dup_samples)} duplicate sample IDs detected "
                f"across batches:", logfile)
            log(dup_samples[["sample_col_name", "detection_rate",
                              "pass_qc"]].to_string(), logfile)
            log("Keeping the instance with higher detection rate.", logfile)
            qc_all = (qc_all
                      .sort_values("detection_rate", ascending=False)
                      .drop_duplicates(subset="sample_col_name", keep="first"))

        passing_samples = set(
            qc_all[qc_all["pass_qc"]]["sample_col_name"]
        )
        log(f"Samples passing QC: {len(passing_samples)}", logfile)

        # ── Load beta matrices batch by batch ─────────────────────────────────
        # Strategy: find probe intersection first (one pass),
        # then load and align (second pass) — avoids holding all in RAM
        log("Pass 1: collecting probe sets per batch...", logfile)
        common_probes_set = None
        empty_batches = []
        n_nonempty = 0
        for bf in beta_files:
            log(f"  Running {bf}", logfile)
            batch_id = os.path.basename(bf).replace(".betas.csv.gz", "")
            log(f"  batch id {batch_id}", logfile)
            try:
                log(f"  Loading csv for {bf}", logfile)
                probes = set(
                    pd.read_csv(bf, index_col=0, usecols=[0]).index
                )
            except Exception as e:
                log(f"  WARNING: could not read {batch_id}: {e}", logfile)
                empty_batches.append(bf)
                continue

            if len(probes) == 0:
                log(f"  SKIPPING {batch_id}: empty batch (no samples passed QC)",
                    logfile)
                empty_batches.append(bf)
                continue

            # Running intersection instead of appending every batch's full
            # set to a list -- same final result (set intersection is
            # associative/commutative), but peak memory is now two sets at
            # a time instead of one per batch (252 x ~700K-900K strings was
            # the likely cause of the climbing memory usage).
            common_probes_set = (
                probes if common_probes_set is None
                else common_probes_set & probes
            )
            n_nonempty += 1
            log(f"  {batch_id}: {len(probes):,} probes "
                f"(running intersection: {len(common_probes_set):,})", logfile)

        log(f"Non-empty batches: {n_nonempty} / {len(beta_files)}",
            logfile)
        if empty_batches:
            log(f"Skipped {len(empty_batches)} empty batches:", logfile)
            for bf in empty_batches:
                log(f"  {os.path.basename(bf)}", logfile)

        if common_probes_set is None:
            log("ERROR: all batches are empty — cannot proceed", logfile)
            sys.exit(1)

        common_probes = sorted(common_probes_set)
        log(f"Common probes across all non-empty batches: {len(common_probes):,}",
            logfile)

        # ── Pass 2: incremental memmap-based concatenation ───────────────────
        # pd.concat holds all dataframes in RAM simultaneously — too expensive
        # for 250 batches x 543K probes. Instead we pre-allocate a numpy
        # memory-mapped array on disk and write each batch into it column by
        # column. Peak RAM = one batch at a time (~50 MB) regardless of total.
        log("Pass 2: building sample list and allocating memmap...", logfile)

        non_empty_files = [bf for bf in beta_files if bf not in empty_batches]

        # Sub-pass 2a: collect all passing sample names in order
        all_samples = []
        valid_files = []
        for bf in non_empty_files:
            batch_id = os.path.basename(bf).replace(".betas.csv.gz", "")
            try:
                cols = pd.read_csv(bf, index_col=0, nrows=0).columns.tolist()
            except Exception as e:
                log(f"  WARNING: could not read header of {batch_id}: {e}",
                    logfile)
                continue
            keep_cols = [c for c in cols if c in passing_samples]
            if not keep_cols:
                log(f"  SKIPPING {batch_id}: no passing samples", logfile)
                continue
            all_samples.extend(keep_cols)
            valid_files.append((bf, keep_cols))

        # Remove duplicates keeping first occurrence (higher detection rate
        # already handled in QC dedup step above)
        seen = set()
        unique_samples = []
        for s in all_samples:
            if s not in seen:
                unique_samples.append(s)
                seen.add(s)

        n_probes  = len(common_probes)
        n_samples = len(unique_samples)
        log(f"Total samples to write: {n_samples}", logfile)
        log(f"Total probes:           {n_probes:,}", logfile)

        if n_samples == 0 or n_probes == 0:
            log("ERROR: no valid samples or probes — cannot proceed", logfile)
            sys.exit(1)

        # Allocate a plain in-memory array. The disk-backed memmap approach
        # was designed to avoid RAM cost for very large matrices, but this
        # one is ~4.2GB -- comfortably within RAM on this VM. A memmap here
        # was actively counterproductive: writing whole *columns*
        # (mm[:, col_indices] = ...) into a row-major disk-backed array
        # means every column write touches scattered, non-contiguous pages
        # across the entire file, which is a pathologically slow I/O
        # pattern for concurrent threaded writes. An in-memory array has no
        # such penalty.
        log(f"Allocating in-memory array "
            f"({n_probes:,} x {n_samples} float32 = "
            f"{n_probes * n_samples * 4 / 1e9:.2f} GB)...", logfile)

        mm = np.full((n_probes, n_samples), np.nan, dtype="float32")

        # Build sample index for fast column placement
        sample_idx = {s: i for i, s in enumerate(unique_samples)}
        probe_idx  = {p: i for i, p in enumerate(common_probes)}

        # Sub-pass 2b: load each batch and write into memmap in parallel
        # Each batch writes to non-overlapping columns so parallel writes are safe.
        # Vectorized column writes replace the per-column Python loop.
        log("Pass 2: writing batches into memmap (parallel)...", logfile)

        def write_batch(args_tuple):
            bf, keep_cols = args_tuple
            batch_id = os.path.basename(bf).replace(".betas.csv.gz", "")
            try:
                df = pd.read_csv(
                    bf,
                    index_col = 0,
                    usecols   = ["probe_id"] + keep_cols
                )
            except Exception as e:
                return batch_id, 0, f"read error: {e}"

            # Align probes in one operation
            df = df.reindex(common_probes)

            # Find column indices for all samples in this batch at once
            valid_cols = [c for c in keep_cols
                          if c in sample_idx and c in df.columns]
            if not valid_cols:
                return batch_id, 0, "no valid columns after filter"

            col_indices = [sample_idx[c] for c in valid_cols]

            # Vectorized write: all columns for this batch in one numpy op
            # mm is a memmap — writes go directly to disk-backed memory
            mm[:, col_indices] = df[valid_cols].values.astype("float32")

            return batch_id, len(valid_cols), None

        n_workers = min(4, len(valid_files))  # limit to 4 to avoid I/O thrash
        n_done    = 0
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(write_batch, item): item
                for item in valid_files
            }
            for future in as_completed(futures):
                batch_id, n_written, err = future.result()
                n_done += 1
                if err:
                    log(f"  WARNING {batch_id}: {err}", logfile)
                else:
                    log(f"  [{n_done}/{len(valid_files)}] "
                        f"{batch_id}: wrote {n_written} samples", logfile)

        # Array is already fully populated in-memory -- no flush needed.
        log("All batches written to in-memory array.", logfile)

        # Build combined DataFrame from the in-memory array for downstream steps
        log("Building final DataFrame from array...", logfile)
        beta_combined = pd.DataFrame(
            mm,
            index   = common_probes,
            columns = unique_samples,
            dtype   = "float32"
        )
        del mm   # release the array once copied into the DataFrame
        log(f"Combined matrix: {beta_combined.shape[0]:,} probes x "
            f"{beta_combined.shape[1]} samples", logfile)

        # ── Final probe-level NA filter across all samples ────────────────────
        na_rate = beta_combined.isna().mean(axis=1)
        n_before = len(beta_combined)
        beta_combined = beta_combined[
            na_rate <= (1 - args.min_probe_detection)
        ]
        log(f"After final NA filter: {n_before:,} -> "
            f"{len(beta_combined):,} probes", logfile)

        # ── Check for unexpected duplicate columns ────────────────────────────
        dup_cols = beta_combined.columns[
            beta_combined.columns.duplicated()
        ].tolist()
        if dup_cols:
            log(f"WARNING: duplicate columns after merge: {dup_cols}",
                logfile)
            beta_combined = beta_combined.loc[
                :, ~beta_combined.columns.duplicated(keep="first")
            ]

        # ── Report any columns with unrecognized visit codes ─────────────────
        bad_cols = [c for c in beta_combined.columns
                    if "_T" in c and not c.rsplit("_T", 1)[-1].isdigit()]
        if bad_cols:
            log(f"WARNING: {len(bad_cols)} columns have non-numeric visit "
                f"suffixes and will be dropped:", logfile)
            for c in bad_cols[:20]:
                log(f"  {c}", logfile)
            beta_combined = beta_combined.drop(columns=bad_cols)
            log(f"Remaining samples after drop: {beta_combined.shape[1]}",
                logfile)

        # ── Sort columns by SubjectID then visit number ───────────────────────
        def sort_key(col):
            parts = col.rsplit("_T", 1)
            subject = parts[0]
            try:
                visit = int(parts[1]) if len(parts) > 1 else 0
            except ValueError:
                visit = 999   # put unrecognized visits at end
            return (subject, visit)

        beta_combined = beta_combined[
            sorted(beta_combined.columns, key=sort_key)
        ]

        # ── Save outputs ──────────────────────────────────────────────────────
        # Write CSV first, then explicitly delete the DataFrame before writing
        # parquet. pa.Table.from_pandas() creates a full in-memory copy which
        # doubles peak RAM. By deleting beta_combined first and re-reading the
        # memmap directly, only one copy is ever in memory at a time.

        # ── Save CSV ──────────────────────────────────────────────────────────
        log(f"Saving combined matrix to {args.out_matrix}...", logfile)
        beta_combined.to_csv(args.out_matrix, compression="gzip")
        log("CSV saved.", logfile)

        # ── Save QC summary and probe list while DataFrame is still in scope ──
        log(f"Saving QC summary to {args.out_qc}...", logfile)
        qc_all.to_csv(args.out_qc, index=False)

        log(f"Saving probe list to {args.out_probes}...", logfile)
        with open(args.out_probes, "w") as fh:
            fh.write("\n".join(beta_combined.index.tolist()) + "\n")

        # ── Final summary while DataFrame is still in scope ───────────────────
        log("=== Final Summary ===", logfile)
        log(f"  Probes retained:   {len(beta_combined):,}", logfile)
        log(f"  Samples retained:  {len(beta_combined.columns)}", logfile)
        log(f"  Matrix size (est): "
            f"{beta_combined.memory_usage(deep=True).sum() / 1e9:.2f} GB",
            logfile)

        visits = pd.Series(beta_combined.columns).str.extract(r"_T(\d+)$")[0]
        log("  Samples per visit:", logfile)
        for v, n in visits.value_counts().sort_index().items():
            log(f"    Visit {v}: {n} samples", logfile)

        # (No memmap temp file to clean up anymore -- array was in-memory only.)


if __name__ == "__main__":
    main()
