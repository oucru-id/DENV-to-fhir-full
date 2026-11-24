#!/usr/bin/env python3
"""Classify Dengue virus serotype from minimap2 PAF alignments against the tagged panel."""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict


def load_panel(panel_path):
    by_acc = {}
    by_serotype = defaultdict(list)
    with open(panel_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            acc = (row.get("accession") or "").strip()
            if not acc:
                continue
            entry = {
                "accession": acc,
                "fasta": row["fasta"].strip(),
                "serotype": row["serotype"].strip(),
                "sylvatic": row.get("sylvatic", "false").strip().lower() == "true",
                "nextclade_dataset": row["nextclade_dataset"].strip(),
                "description": row.get("description", "").strip(),
            }
            by_acc[acc] = entry
            by_serotype[entry["serotype"]].append(entry)
    return by_acc, by_serotype


def parse_target_name(tname, by_acc):
    """Extract (serotype, accession, sylvatic) from a tagged target like 'DENV-1|NC_001477.1'."""
    if "|" in tname:
        _tag, acc = tname.split("|", 1)
    else:
        acc = tname
    if acc in by_acc:
        entry = by_acc[acc]
        return entry["serotype"], entry["accession"], entry["sylvatic"]
    # Fallback match by accession prefix (version-less)
    base = acc.split(".")[0]
    for known_acc, entry in by_acc.items():
        if known_acc.split(".")[0] == base:
            return entry["serotype"], entry["accession"], entry["sylvatic"]
    return None, acc, False


def classify_paf(paf_path, panel_path, sample_id, platform, total_reads,
                 min_reads=10, coinfection_min_fraction=0.10):
    by_acc, _by_sero = load_panel(panel_path)
    min_alen = 200 if platform == "fasta" else 50
    min_ident = 0.75 if platform == "fasta" else 0.70

    best_per_query = {}
    if paf_path and os.path.exists(paf_path) and os.path.getsize(paf_path) > 0:
        with open(paf_path, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 12:
                    continue
                qname = parts[0]
                tname = parts[5]
                try:
                    nmatch = int(parts[9])
                    alen = int(parts[10])
                    mapq = int(parts[11])
                except ValueError:
                    continue
                if alen < min_alen or alen <= 0:
                    continue
                ident = nmatch / alen
                if ident < min_ident:
                    continue
                serotype, acc, sylvatic = parse_target_name(tname, by_acc)
                if not serotype:
                    continue
                prev = best_per_query.get(qname)
                if prev is None or (nmatch, alen, mapq) > (prev["nmatch"], prev["alen"], prev["mapq"]):
                    best_per_query[qname] = {
                        "serotype": serotype,
                        "accession": acc,
                        "sylvatic": sylvatic,
                        "nmatch": nmatch,
                        "alen": alen,
                        "identity": ident,
                        "mapq": mapq,
                    }

    sero_counts = defaultdict(int)
    sero_idents = defaultdict(list)
    acc_counts = defaultdict(int)
    for hit in best_per_query.values():
        sero_counts[hit["serotype"]] += 1
        sero_idents[hit["serotype"]].append(hit["identity"])
        acc_counts[hit["accession"]] += 1

    mapped_reads = len(best_per_query)
    ranked_seros = sorted(sero_counts.items(), key=lambda kv: (-kv[1], kv[0]))

    breakdown = {}
    for sero, count in ranked_seros:
        breakdown[sero] = {
            "reads": count,
            "fraction": round(count / mapped_reads, 4) if mapped_reads else 0.0,
            "mean_identity": round(sum(sero_idents[sero]) / len(sero_idents[sero]), 4),
        }

    result = {
        "sample_id": sample_id,
        "platform": platform,
        "method": "minimap2_panel",
        "total_input_reads": total_reads,
        "mapped_denv_reads": mapped_reads,
        "serotype": "unknown",
        "confidence": "none",
        "sylvatic": False,
        "selected_accession": None,
        "selected_reference": "UNKNOWN",
        "nextclade_dataset": "unknown",
        "primary_fraction": 0.0,
        "primary_mean_identity": 0.0,
        "coinfection_suspected": False,
        "secondary_serotype": None,
        "secondary_reads": 0,
        "secondary_fraction": 0.0,
        "serotype_breakdown": breakdown,
        "accession_counts": dict(sorted(acc_counts.items(), key=lambda kv: (-kv[1], kv[0]))),
    }

    if not ranked_seros or ranked_seros[0][1] < min_reads:
        result["note"] = (
            f"insufficient_denv_reads ({mapped_reads} mapped < min_reads {min_reads})"
            if mapped_reads > 0 else "no_denv_reads_mapped"
        )
        return result

    top_sero, top_count = ranked_seros[0]
    top_frac = top_count / mapped_reads
    top_ident = sum(sero_idents[top_sero]) / len(sero_idents[top_sero])

    # Choose the best reference accession belonging to the winning serotype
    sero_accs = [(acc, cnt) for acc, cnt in acc_counts.items() if by_acc[acc]["serotype"] == top_sero]
    sero_accs.sort(key=lambda kv: (-kv[1], by_acc[kv[0]]["sylvatic"], kv[0]))
    best_acc = sero_accs[0][0]
    ref_entry = by_acc[best_acc]

    # Co-infection check
    coinfection = False
    sec_sero = None
    sec_count = 0
    sec_frac = 0.0
    if len(ranked_seros) > 1 and platform != "fasta":
        sec_sero, sec_count = ranked_seros[1]
        sec_frac = sec_count / mapped_reads
        if sec_count >= min_reads and sec_frac >= coinfection_min_fraction:
            coinfection = True

    if platform == "fasta":
        confidence = "high" if top_ident >= 0.88 else ("medium" if top_ident >= 0.80 else "low")
    else:
        if top_count >= 100 and top_frac >= 0.90 and top_ident >= 0.85:
            confidence = "high"
        elif top_count >= 25 and top_frac >= 0.75 and top_ident >= 0.80:
            confidence = "medium"
        else:
            confidence = "low"
        if coinfection and confidence == "high":
            confidence = "medium"

    result.update({
        "serotype": top_sero,
        "confidence": confidence,
        "sylvatic": ref_entry["sylvatic"],
        "selected_accession": ref_entry["accession"],
        "selected_reference": ref_entry["fasta"],
        "nextclade_dataset": ref_entry["nextclade_dataset"],
        "primary_fraction": round(top_frac, 4),
        "primary_mean_identity": round(top_ident, 4),
        "coinfection_suspected": coinfection,
        "secondary_serotype": sec_sero if coinfection else None,
        "secondary_reads": sec_count if coinfection else 0,
        "secondary_fraction": round(sec_frac, 4) if coinfection else 0.0,
    })
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", required=True, help="minimap2 PAF alignment file")
    parser.add_argument("--panel", required=True, help="reference panel TSV")
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--platform", default="illumina", choices=["illumina", "ont", "fasta"])
    parser.add_argument("--total_reads", type=int, default=0)
    parser.add_argument("--min_reads", type=int, default=10)
    parser.add_argument("--coinfection_min_fraction", type=float, default=0.10)
    parser.add_argument("--output", required=True, help="output serotype JSON path")
    parser.add_argument("--env_output", default=None, help="optional shell env assignment file")
    args = parser.parse_args(argv)

    result = classify_paf(
        paf_path=args.paf,
        panel_path=args.panel,
        sample_id=args.sample_id,
        platform=args.platform,
        total_reads=args.total_reads,
        min_reads=args.min_reads,
        coinfection_min_fraction=args.coinfection_min_fraction,
    )

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)

    if args.env_output:
        with open(args.env_output, "w", encoding="utf-8") as handle:
            handle.write(f'DENV_SEROTYPE="{result["serotype"]}"\n')
            handle.write(f'DENV_REF_FILE="{result["selected_reference"]}"\n')
            handle.write(f'DENV_NEXTCLADE_DATASET="{result["nextclade_dataset"]}"\n')

    print(
        f"[{args.sample_id}] serotype={result['serotype']} "
        f"ref={result['selected_reference']} mapped={result['mapped_denv_reads']} "
        f"confidence={result['confidence']} coinfection={result['coinfection_suspected']}"
    )


if __name__ == "__main__":
    main()
