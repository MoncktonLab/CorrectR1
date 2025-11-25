#!/usr/bin/env python3
"""
CorrectR1_v1.12.py
Parallel R1 correction across multiple samples in a directory.

Outputs organized as:
  R1corrected_<params>_reads/
      sample1_correctedR1_R1.fastq.gz
      sample1_correctedR1_R2.fastq.gz
      ...
  R1corrected_<params>_figures/
      sample1_correction_counts.png
      sample1_correction_positions.png
      sample1_overlap_heatmap.png
      ...
  combined_summary.tsv
  combined_run.log
"""

import argparse
import gzip
import os
import glob
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed, ProcessPoolExecutor
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd
import time
from datetime import datetime

# --- Helpers ---
def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def open_maybe_gz(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def rc(seq):
    return str(Seq(seq).reverse_complement())

def find_best_local_overlap(r1_seq, r2_rc_seq, aligner, min_ov=20):
    if not r1_seq or not r2_rc_seq:
        return None
    alignments = aligner.align(r1_seq, r2_rc_seq)
    if not alignments:
        return None
    best = alignments[0]
    r1_start, r1_end = best.aligned[0][0]
    r2_start, r2_end = best.aligned[1][0]
    overlap_len = r1_end - r1_start
    if overlap_len < min_ov:
        return None
    return {'r1_start': r1_start, 'r1_end': r1_end, 'r2_start': r2_start, 'r2_end': r2_end}

def correct_pair(r1_rec, r2_rec, aligner, r2_rc_trim_start=0, r2_rc_trim_end=0, min_ov=20):
    r1_seq = str(r1_rec.seq)
    r1_q = r1_rec.letter_annotations.get("phred_quality", [])
    r2_rc_seq = rc(r2_rec.seq)
    r2_q = list(reversed(r2_rec.letter_annotations.get("phred_quality", [])))

    # Handle trimming
    if r2_rc_trim_end == 0:
        if r2_rc_trim_start > 0:
            r2_rc_seq = r2_rc_seq[-r2_rc_trim_start:]
            r2_q = r2_q[-r2_rc_trim_start:]
    else:
        if r2_rc_trim_start > 0:
            r2_rc_seq = r2_rc_seq[-r2_rc_trim_start:-r2_rc_trim_end]
            r2_q = r2_q[-r2_rc_trim_start:-r2_rc_trim_end]
        else:
            r2_rc_seq = r2_rc_seq[:-r2_rc_trim_end]
            r2_q = r2_q[:-r2_rc_trim_end]

    if len(r1_seq) == 0 or len(r2_rc_seq) == 0:
        return r1_rec, 0, [], None, None

    ov = find_best_local_overlap(r1_seq, r2_rc_seq, aligner, min_ov)
    if ov is None:
        return r1_rec, 0, [], None, None

    r1_start, r1_end, r2_start, r2_end = ov['r1_start'], ov['r1_end'], ov['r2_start'], ov['r2_end']
    overlap_len = r1_end - r1_start

    r1_seq_list, r1_qual_list = list(r1_seq), list(r1_q)
    r2_overlap_seq = r2_rc_seq[r2_start:r2_end]
    r2_overlap_qual = r2_q[r2_start:r2_end]

    bases_corrected = 0
    corrected_positions = []

    for i in range(overlap_len):
        idx = r1_start + i
        if idx >= len(r1_seq_list): break
        b1, q1 = r1_seq_list[idx], r1_qual_list[idx] if idx < len(r1_qual_list) else 0
        b2, q2 = r2_overlap_seq[i], r2_overlap_qual[i] if i < len(r2_overlap_qual) else 0
        if b1 != b2 and q2 > q1:
            r1_seq_list[idx], r1_qual_list[idx] = b2, q2
            bases_corrected += 1
            corrected_positions.append(idx)

    r1_rec.seq = Seq(''.join(r1_seq_list))
    r1_rec.letter_annotations["phred_quality"] = r1_qual_list
    return r1_rec, bases_corrected, corrected_positions, r1_start, r1_end


# --- Main per-sample function ---
def process_sample(r1_path, r2_path, sample_name, reads_dir, figs_dir,
                   skip_n, r2_rc_trim_start, r2_rc_trim_end,
                   match_score, mismatch_score, open_gap_score, extend_gap_score, motif):
    start_time = time.time()
    log_lines = []
    log_lines.append(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Starting sample correction for {r1_path}")

    # --- Alignment setup ---
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    # --- Output paths ---
    r1_out = os.path.join(reads_dir, f"{sample_name}_correctedR1_R1.fastq.gz")
    r2_out = os.path.join(reads_dir, f"{sample_name}_correctedR1_R2.fastq.gz")

    # --- Read input ---
    r1_handle = open_maybe_gz(r1_path, 'rt')
    r2_handle = open_maybe_gz(r2_path, 'rt')
    r1_iter = islice(SeqIO.parse(r1_handle, "fastq"), skip_n, None)
    r2_iter = islice(SeqIO.parse(r2_handle, "fastq"), skip_n, None)

    corrected_counts, corrected_positions_all, overlap_starts, overlap_ends = [], [], [], []
    total, max_r1_len, pass_reads_before = 0, 0, 0

    # --- Correction loop ---
    with gzip.open(r1_out, "wt") as r1_out_handle, gzip.open(r2_out, "wt") as r2_out_handle:
        for r1, r2 in zip(r1_iter, r2_iter):
            if motif in str(r1.seq):
                pass_reads_before += 1
            r1_rec, bases_corrected, corrected_positions, r1_start, r1_end = correct_pair(
                r1, r2, aligner,
                r2_rc_trim_start=r2_rc_trim_start,
                r2_rc_trim_end=r2_rc_trim_end,
                min_ov=20
            )
            SeqIO.write(r1_rec, r1_out_handle, "fastq")
            SeqIO.write(r2, r2_out_handle, "fastq")
            corrected_counts.append(bases_corrected)
            corrected_positions_all.extend(corrected_positions)
            if r1_start is not None and r1_end is not None:
                overlap_starts.append(r1_start)
                overlap_ends.append(r1_end)
            max_r1_len = max(max_r1_len, len(r1_rec.seq))
            total += 1

    r1_handle.close()
    r2_handle.close()

    # --- Count motif occurrences ---
    pass_reads_after = sum(1 for rec in SeqIO.parse(gzip.open(r1_out, "rt"), "fastq") if motif in str(rec.seq))

    # --- Plot generation (isolated per process, thread-safe) ---
    if corrected_counts:
        plt.figure(figsize=(8,5))
        plt.bar(*zip(*Counter(corrected_counts).items()))
        plt.xlabel("Bases corrected per read")
        plt.ylabel("Number of reads")
        plt.title(f"{sample_name} corrections per read")
        plt.tight_layout()
        plt.savefig(os.path.join(figs_dir, f"{sample_name}_correction_counts.png"), dpi=150)
        plt.close()

    if corrected_positions_all:
        plt.figure(figsize=(10,5))
        plt.hist(corrected_positions_all, bins=range(0, max(corrected_positions_all)+5, 5))
        plt.xlabel("R1 position")
        plt.ylabel("Corrections")
        plt.title(f"{sample_name} positional correction distribution")
        plt.tight_layout()
        plt.savefig(os.path.join(figs_dir, f"{sample_name}_correction_positions.png"), dpi=150)
        plt.close()

    if overlap_starts and overlap_ends:
        coverage = np.zeros(max_r1_len, dtype=int)
        for s, e in zip(overlap_starts, overlap_ends):
            if s < e and e <= max_r1_len:
                coverage[s:e] += 1
        coverage_fraction = coverage / len(overlap_starts)
        plt.figure(figsize=(12,2))
        plt.imshow([coverage_fraction], aspect="auto", cmap="viridis", extent=[0, max_r1_len, 0, 1])
        plt.colorbar(label="Fraction of reads with overlap")
        plt.xlabel("R1 position")
        plt.yticks([])
        plt.title(f"{sample_name} overlap coverage")
        plt.tight_layout()
        plt.savefig(os.path.join(figs_dir, f"{sample_name}_overlap_heatmap.png"), dpi=150)
        plt.close()

    # --- Prepare summary row ---
    bases_counter = Counter(corrected_counts)
    max_correction = 30
    bases_columns = [f"{i}_bases_corrected" for i in range(max_correction + 1)]
    bases_values = [bases_counter.get(i, 0) for i in range(max_correction + 1)]

    pos_counter = Counter(corrected_positions_all)
    max_pos = max(pos_counter.keys()) if pos_counter else 0
    pos_columns = [f"pos_{i}" for i in range(max_pos + 1)]
    pos_values = [pos_counter.get(i, 0) for i in range(max_pos + 1)]

    end_time = time.time()
    elapsed = end_time - start_time  

    # --- Summary row dictionary ---
    row = {
        "sample": sample_name,
        "total_pairs": total,
        "pass_reads_before": pass_reads_before,
        "pass_reads_after": pass_reads_after,
        "match_score": match_score,
        "mismatch_score": mismatch_score,
        "gapopen": open_gap_score,
        "gapextend": extend_gap_score,
        "trimstart": r2_rc_trim_start,
        "trimend": r2_rc_trim_end,
        "mean_corrections": np.mean(corrected_counts) if corrected_counts else 0,
        "timetaken": round(elapsed, 2)
    }
    row.update(dict(zip(bases_columns, bases_values)))
    row.update(dict(zip(pos_columns, pos_values)))

    # --- Add summary info to log ---
    if corrected_counts:
        log_lines.append(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Summary of bases corrected per read: "
                         f"mean={np.mean(corrected_counts):.2f}, max={max(corrected_counts)}, "
                         f"median={np.median(corrected_counts):.2f}, total_reads={len(corrected_counts)}")
        corr_dist = ", ".join([f"{k}:{v}" for k, v in sorted(Counter(corrected_counts).items())])
        log_lines.append(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Correction distribution (bases_corrected:count): {corr_dist}")

    if corrected_positions_all:
        pos_min, pos_max = min(corrected_positions_all), max(corrected_positions_all)
        log_lines.append(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Corrected positions range: {pos_min}-{pos_max}, total_positions={len(corrected_positions_all)}")

    log_lines.append(f"Motif occurrences before correction: {pass_reads_before}")
    log_lines.append(f"Motif occurrences after correction: {pass_reads_after}")
    
    log_lines.append(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Finished sample {sample_name} in {elapsed:.2f} seconds.")

    return "\n".join(log_lines), row


# --- Entry point ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch R1 correction using R2 overlap")
    parser.add_argument('-d', '--dir', required=True, help='Directory containing paired FASTQ files')
    parser.add_argument('-t', type=int, default=1, help='Number of threads (samples in parallel)')
    parser.add_argument('-s', type=int, default=0, help='Reads to skip from start')
    parser.add_argument('--trimstart', type=int, default=0)
    parser.add_argument('--trimend', type=int, default=0)
    parser.add_argument('-m', type=int, default=2)
    parser.add_argument('-mm', type=int, default=-1)
    parser.add_argument('-go', type=int, default=-5)
    parser.add_argument('-ge', type=int, default=-1)
    parser.add_argument('--motif', type=str, default="CAACAGCCGCCA")
    args = parser.parse_args()

    # Parameter-specific subfolders
    cwd = os.getcwd()
    batch = os.path.basename(args.dir)
    param_suffix = f"m{args.m}_mm{args.mm}_go{args.go}_ge{args.ge}_ts{args.trimstart}_te{args.trimend}"
    reads_dir = os.path.join(cwd, f"{batch}_R1corrected_{param_suffix}_reads")
    figs_dir = os.path.join(cwd, f"{batch}_R1corrected_{param_suffix}_figures")
    os.makedirs(reads_dir, exist_ok=True)
    os.makedirs(figs_dir, exist_ok=True)
    
    # Find paired samples
    r1_files = sorted(glob.glob(os.path.join(args.dir, "*_R1*.fastq*")))
    jobs = []
    for r1 in r1_files:
        sample_name = os.path.basename(r1).split("_R1")[0]
        r2 = r1.replace("_R1", "_R2")
        if not os.path.exists(r2):
            print(f"Missing R2 for {sample_name}")
            continue
        jobs.append((r1, r2, sample_name))

    print(f"Found {len(jobs)} samples. Running up to {args.t} in parallel...\n")

    summary_rows, logs = [], []
    
    with ProcessPoolExecutor(max_workers=args.t) as ex:
        future_to_sample = {
            ex.submit(
                process_sample,
                r1, r2, sample, reads_dir, figs_dir,
                args.s, args.trimstart, args.trimend,
                args.m, args.mm, args.go, args.ge, args.motif
            ): sample
            for (r1, r2, sample) in jobs
        }
    
        for f in as_completed(future_to_sample):
            sname = future_to_sample[f]
            try:
                sample_log, summary_row = f.result()  
                logs.append(sample_log)
                summary_rows.append(summary_row)
                print(f"{sname} done.")
            except Exception as e:
                logs.append(f"[ERROR] {sname}: {e}")
                print(f"{sname} failed: {e}")

    # Write combined log
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pythonscriptname = os.path.basename(__file__)
    log_path = os.path.join(cwd, f"{batch}_MiSeq_CorrectR1_combined_run_{timestamp}.log")
    with open(log_path, "w") as outfile:
        outfile.write(f"# --- Samples in {args.dir} processed using {pythonscriptname} at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ---\n\n")
        for sample_log in logs:
            outfile.write(f"# --- Sample log at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ---\n")
            outfile.write(sample_log)
            outfile.write("\n\n")

    # --- Write combined summary TSV ---
    summary_path = os.path.join(cwd, f"{batch}_MiSeq_CorrectR1_combined_summary_{timestamp}.tsv")
    if summary_rows:
        df = pd.DataFrame(summary_rows)
        df.to_csv(summary_path, sep='\t', index=False)
    else:
        print("Warning: No summary rows collected. TSV not written.")

    print(f"\nSummary: {summary_path}")
    print(f"Log: {log_path}")
    print("All samples complete.\n")



## Example usage  
#python correct_R1_v1.12.py -d "/home/mw304m/CorrectR1/HTT_Sperm/PanelC" -t 8 --trimstart 50 --trimend 30 


# python correct_R1_v1.12.py -d "DATA-00001311 - MiSeq data consiting of amplicon sequencing data of the HTT exon 1 repeat locus/Danish_HTT_TargetedMiSeq_Reads" -t 8 --trimstart 50 --trimend 30 