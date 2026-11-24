#!/usr/bin/env python3
"""Generate human-readable per-sample Dengue genomic summary report."""

import argparse
import json
from collections import defaultdict
from datetime import datetime, timezone

from summarize_coverage import read_first_fasta, sequence_stats


def format_report(sample_id, consensus_path, serotype_data, genotype_data, coverage_data):
    _hdr, seq = read_first_fasta(consensus_path)
    seq_stats = genotype_data.get("sequence_stats") or sequence_stats(seq)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines = [
        "=" * 80,
        "DENGUE VIRUS GENOMIC ANALYSIS REPORT (DENVtoFHIR 2.0)",
        f"Sample ID: {sample_id}",
        f"Generated: {timestamp}",
        "=" * 80,
        "",
        "SEQUENCE & COVERAGE QUALITY:",
        "-" * 50,
        f"Consensus length               : {seq_stats['length']:,} bp",
        f"Consensus completeness (non-N) : {seq_stats['completeness'] * 100:.2f}%",
        f"Ambiguous bases (N)            : {seq_stats['n_bases']:,} ({seq_stats['n_percentage']:.2f}%)",
        f"IUPAC non-ACGTN bases          : {seq_stats['non_acgtn_bases']:,}",
        f"GC content                     : {seq_stats['gc_percentage']:.2f}%",
    ]

    if coverage_data.get("assessed"):
        min_dp = coverage_data.get("min_depth_threshold", 10)
        lines.extend([
            f"Reference genome               : {coverage_data.get('reference')}",
            f"Genome breadth (>=1x)          : {(coverage_data.get('breadth_1x') or 0) * 100:.2f}%",
            f"Genome breadth (>={min_dp}x, masked) : {(coverage_data.get('breadth_min_depth') or 0) * 100:.2f}%",
            f"Genome breadth (>=30x)         : {(coverage_data.get('breadth_30x') or 0) * 100:.2f}%",
            f"Mean sequencing depth          : {coverage_data.get('mean_depth', 0):.1f}x",
            f"Median sequencing depth        : {coverage_data.get('median_depth', 0):.1f}x",
            f"Low-depth masked bases (<{min_dp}x)  : {coverage_data.get('masked_bases', 0):,}",
        ])
    else:
        lines.append(f"Coverage assessment            : {coverage_data.get('status', 'not_assessed')}")

    lines.extend([
        "",
        "SEROTYPE CLASSIFICATION:",
        "-" * 50,
        f"Serotype                       : {serotype_data.get('serotype', 'unknown')}",
        f"Confidence                     : {serotype_data.get('confidence', 'none')}",
        f"Selected reference             : {serotype_data.get('selected_accession') or 'none'}",
        f"Sylvatic reference match       : {'YES' if serotype_data.get('sylvatic') else 'no'}",
        f"Mapped DENV reads              : {serotype_data.get('mapped_denv_reads', 0):,} / {serotype_data.get('total_input_reads', 0):,}",
        f"Primary serotype fraction      : {(serotype_data.get('primary_fraction') or 0) * 100:.2f}%",
    ])
    if serotype_data.get("coinfection_suspected"):
        lines.append(
            f"WARNING - Co-infection flagged : {serotype_data.get('secondary_serotype')} "
            f"({serotype_data.get('secondary_reads')} reads, "
            f"{(serotype_data.get('secondary_fraction') or 0) * 100:.2f}%)"
        )

    nc = genotype_data.get("nextclade", {})
    lines.extend([
        "",
        "GENOTYPE & LINEAGE CLASSIFICATION (Nextclade):",
        "-" * 50,
        f"QC status                      : {genotype_data.get('qc_status', 'UNKNOWN')}",
        f"Genotype                       : {genotype_data.get('genotype', 'Unknown')}",
        f"Major lineage                  : {genotype_data.get('major_lineage', 'Unknown')}",
        f"Minor lineage                  : {genotype_data.get('minor_lineage', 'Unknown')}",
        f"Clade code                     : {genotype_data.get('clade', 'Unknown')}",
        f"Confidence                     : {genotype_data.get('confidence', 'none')}",
        f"Nextclade QC status (score)    : {nc.get('qc_status', 'not_run')} ({nc.get('qc_score')})",
    ])
    if genotype_data.get("qc_reasons"):
        lines.append(f"QC flags / reasons             : {', '.join(genotype_data['qc_reasons'])}")
        if genotype_data.get("raw_clade") not in (None, "", "unknown"):
            lines.append(f"Unfiltered Nextclade clade     : {genotype_data['raw_clade']}")
    if nc.get("frameshifts"):
        lines.append(f"Frameshifts                    : {', '.join(nc['frameshifts'])}")

    aa_subs = genotype_data.get("aa_substitutions", [])
    aa_dels = genotype_data.get("aa_deletions", [])
    aa_ins = genotype_data.get("aa_insertions", [])

    lines.extend([
        "",
        f"DETECTED AMINO ACID CHANGES (substitutions={len(aa_subs)}, deletions={len(aa_dels)}, insertions={len(aa_ins)}):",
        "-" * 50,
    ])

    by_gene = defaultdict(list)
    for m in aa_subs:
        label = f"{m['ref_aa']}{m['codon']}{m['alt_aa']}" if m.get("codon") else m["raw"]
        by_gene[m["gene"]].append(label)
    for m in aa_dels:
        label = f"{m['ref_aa']}{m['codon']}del" if m.get("codon") else m["raw"]
        by_gene[m["gene"]].append(label)
    for m in aa_ins:
        label = f"{m['codon']}ins{m['inserted_aa']}" if m.get("codon") else m["raw"]
        by_gene[m["gene"]].append(label)

    if by_gene:
        for gene in sorted(by_gene):
            lines.append(f"  {gene}: {', '.join(by_gene[gene])}")
    else:
        lines.append("  None detected")

    lines.extend(["", "=" * 80, "END OF DENGUE ANALYSIS REPORT", "=" * 80, ""])
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--serotype_json", required=True)
    parser.add_argument("--genotype_json", required=True)
    parser.add_argument("--coverage_json", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    with open(args.serotype_json, encoding="utf-8") as handle:
        serotype_data = json.load(handle)
    with open(args.genotype_json, encoding="utf-8") as handle:
        genotype_data = json.load(handle)
    with open(args.coverage_json, encoding="utf-8") as handle:
        coverage_data = json.load(handle)

    text = format_report(args.sample_id, args.consensus, serotype_data, genotype_data, coverage_data)
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(text)


if __name__ == "__main__":
    main()
