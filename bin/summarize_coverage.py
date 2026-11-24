#!/usr/bin/env python3
"""Summarize per-base coverage from mosdepth, generate low-depth mask BED, or inspect FASTA inputs."""

import argparse
import gzip
import json
import os
import sys


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def read_first_fasta(path):
    header = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:].split()[0]
            else:
                chunks.append(line.upper())
    return header or "", "".join(chunks)


def sequence_stats(seq):
    length = len(seq)
    if length == 0:
        return {
            "length": 0,
            "non_n_bases": 0,
            "n_bases": 0,
            "non_acgtn_bases": 0,
            "completeness": 0.0,
            "n_percentage": 100.0,
            "gc_percentage": 0.0,
        }
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    n = seq.count("N")
    acgt = a + c + g + t
    non_acgtn = length - acgt - n
    non_n = length - n
    gc = ((g + c) / acgt * 100.0) if acgt > 0 else 0.0
    return {
        "length": length,
        "non_n_bases": non_n,
        "n_bases": n,
        "non_acgtn_bases": non_acgtn,
        "completeness": round(non_n / length, 4),
        "n_percentage": round(n / length * 100.0, 2),
        "gc_percentage": round(gc, 2),
    }


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1][1] = end
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def summarize_mosdepth(sample_id, per_base_bed, reference_fasta, min_depth, output_mask_bed):
    ref_acc, ref_seq = read_first_fasta(reference_fasta)
    ref_len = len(ref_seq)
    if ref_len == 0:
        sys.exit(f"ERROR: empty reference FASTA: {reference_fasta}")

    covered_intervals = []
    depth_runs = []  # (depth, length)
    mask_intervals = []
    bases_ge_1 = 0
    bases_ge_min = 0
    bases_ge_30 = 0
    total_depth_sum = 0
    max_depth = 0
    seen_chrom = ref_acc

    with open_text(per_base_bed) as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = max(0, int(parts[1]))
            end = min(ref_len, int(parts[2]))
            depth = int(float(parts[3]))
            if end <= start:
                continue
            seen_chrom = chrom
            span = end - start
            covered_intervals.append((start, end))
            depth_runs.append((depth, span))
            total_depth_sum += depth * span
            if depth > max_depth:
                max_depth = depth
            if depth >= 1:
                bases_ge_1 += span
            if depth >= min_depth:
                bases_ge_min += span
            if depth >= 30:
                bases_ge_30 += span
            if depth < min_depth:
                mask_intervals.append((start, end))

    # Any reference bases absent from per-base.bed have depth 0 and must be masked
    merged_cov = merge_intervals(covered_intervals)
    cursor = 0
    missing_zero_bases = 0
    for s, e in merged_cov:
        if s > cursor:
            mask_intervals.append((cursor, s))
            missing_zero_bases += s - cursor
        cursor = max(cursor, e)
    if cursor < ref_len:
        mask_intervals.append((cursor, ref_len))
        missing_zero_bases += ref_len - cursor
    if missing_zero_bases > 0:
        depth_runs.append((0, missing_zero_bases))

    merged_mask = merge_intervals(mask_intervals)
    masked_bases = sum(e - s for s, e in merged_mask)

    # Exact weighted median depth across ref_len bases
    depth_runs.sort(key=lambda pair: pair[0])
    min_observed = depth_runs[0][0] if depth_runs else 0
    mid = ref_len / 2.0
    cum = 0
    median_depth = 0.0
    for depth, span in depth_runs:
        cum += span
        if cum >= mid:
            median_depth = float(depth)
            break

    if output_mask_bed:
        with open(output_mask_bed, "w", encoding="utf-8") as out:
            for s, e in merged_mask:
                out.write(f"{seen_chrom}\t{s}\t{e}\n")

    return {
        "sample_id": sample_id,
        "assessed": True,
        "status": "assessed",
        "reference": seen_chrom,
        "reference_length": ref_len,
        "min_depth_threshold": min_depth,
        "mean_depth": round(total_depth_sum / ref_len, 2),
        "median_depth": round(median_depth, 2),
        "min_depth": min_observed,
        "max_depth": max_depth,
        "bases_ge_1x": bases_ge_1,
        "bases_ge_min_depth": bases_ge_min,
        "bases_ge_30x": bases_ge_30,
        "breadth_1x": round(bases_ge_1 / ref_len, 4),
        "breadth_min_depth": round(bases_ge_min / ref_len, 4),
        "breadth_30x": round(bases_ge_30 / ref_len, 4),
        "masked_bases": masked_bases,
        "masked_intervals": len(merged_mask),
    }


def prepare_fasta(sample_id, input_fasta, output_fasta):
    orig_header, seq = read_first_fasta(input_fasta)
    with open(output_fasta, "w", encoding="utf-8") as out:
        out.write(f">{sample_id}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i + 80] + "\n")
    stats = sequence_stats(seq)
    return {
        "sample_id": sample_id,
        "assessed": False,
        "status": "submitted_consensus",
        "original_header": orig_header,
        "reference": None,
        "reference_length": stats["length"],
        "min_depth_threshold": None,
        "mean_depth": None,
        "median_depth": None,
        "breadth_1x": None,
        "breadth_min_depth": None,
        "breadth_30x": None,
        "masked_bases": stats["n_bases"],
        "sequence_stats": stats,
    }


def make_unmapped(sample_id, output_fasta):
    with open(output_fasta, "w", encoding="utf-8") as out:
        out.write(f">{sample_id}\n\n")
    return {
        "sample_id": sample_id,
        "assessed": False,
        "status": "unmapped",
        "reference": None,
        "reference_length": 0,
        "min_depth_threshold": None,
        "mean_depth": 0.0,
        "median_depth": 0.0,
        "breadth_1x": 0.0,
        "breadth_min_depth": 0.0,
        "breadth_30x": 0.0,
        "masked_bases": 0,
        "sequence_stats": sequence_stats(""),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=["mosdepth", "fasta", "unmapped"])
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--per_base_bed", default=None)
    parser.add_argument("--reference", default=None)
    parser.add_argument("--min_depth", type=int, default=10)
    parser.add_argument("--input_fasta", default=None)
    parser.add_argument("--output_fasta", default=None)
    parser.add_argument("--output_mask_bed", default=None)
    parser.add_argument("--output_json", required=True)
    args = parser.parse_args(argv)

    if args.mode == "mosdepth":
        if not args.per_base_bed or not args.reference:
            sys.exit("ERROR: --mode mosdepth requires --per_base_bed and --reference")
        summary = summarize_mosdepth(
            args.sample_id, args.per_base_bed, args.reference,
            args.min_depth, args.output_mask_bed
        )
    elif args.mode == "fasta":
        if not args.input_fasta or not args.output_fasta:
            sys.exit("ERROR: --mode fasta requires --input_fasta and --output_fasta")
        summary = prepare_fasta(args.sample_id, args.input_fasta, args.output_fasta)
    else:
        if not args.output_fasta:
            sys.exit("ERROR: --mode unmapped requires --output_fasta")
        summary = make_unmapped(args.sample_id, args.output_fasta)

    with open(args.output_json, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(
        f"[{args.sample_id}] mode={args.mode} status={summary['status']} "
        f"mean_depth={summary.get('mean_depth')} breadth_min={summary.get('breadth_min_depth')}"
    )


if __name__ == "__main__":
    main()
