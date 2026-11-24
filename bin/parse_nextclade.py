#!/usr/bin/env python3
"""Parse Nextclade CSV output, apply consensus completeness/QC gating, and emit structured genotype JSON."""

import argparse
import csv
import json
import os
import re

from summarize_coverage import read_first_fasta, sequence_stats

# Hill et al. (2024) Dengue lineage nomenclature: e.g. 1V_D.1.1, 2III_C.2, 3III_B.3.2, 4II_B.2, 1I
CLADE_RE = re.compile(r"^([1-4])([IVX]+)(?:[._]([A-Z0-9]+))?(?:[._](.+))?$")
AA_SUB_RE = re.compile(r"^([^:]+):([A-Za-z*])(\d+)([A-Za-z*])$")
AA_DEL_RE = re.compile(r"^([^:]+):([A-Za-z*])(\d+)-$")
AA_INS_RE = re.compile(r"^([^:]+):(\d+):([A-Za-z*]+)$")


def parse_clade_label(clade):
    """Decompose a Nextclade Dengue clade label into genotype, major lineage, and minor lineage."""
    if not clade or clade.lower() in ("unknown", "unassigned", "na"):
        return "Unclassified", None, None
    match = CLADE_RE.match(clade.strip())
    if match:
        geno = f"Genotype {match.group(2)}"
        major = f"Lineage {match.group(3)}" if match.group(3) else None
        minor = f"Sub-lineage {match.group(4)}" if match.group(4) else None
        return geno, major, minor
    if clade.upper().startswith("DENV"):
        return f"{clade} (Basal/Unclassified)", None, None
    return clade, None, None


def split_csv_list(value):
    if not value:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_aa_mutations(row):
    substitutions = []
    for item in split_csv_list(row.get("aaSubstitutions", "")):
        m = AA_SUB_RE.match(item)
        if m:
            substitutions.append({
                "gene": m.group(1),
                "ref_aa": m.group(2),
                "codon": int(m.group(3)),
                "alt_aa": m.group(4),
                "raw": item,
                "hgvs_p": f"p.{m.group(2)}{m.group(3)}{m.group(4)}",
                "type": "substitution",
            })
        elif ":" in item:
            gene, change = item.split(":", 1)
            substitutions.append({
                "gene": gene, "ref_aa": "", "codon": None, "alt_aa": "",
                "raw": item, "hgvs_p": f"p.{change}", "type": "substitution",
            })

    deletions = []
    for item in split_csv_list(row.get("aaDeletions", "")):
        m = AA_DEL_RE.match(item)
        if m:
            deletions.append({
                "gene": m.group(1),
                "ref_aa": m.group(2),
                "codon": int(m.group(3)),
                "alt_aa": "-",
                "raw": item,
                "hgvs_p": f"p.{m.group(2)}{m.group(3)}del",
                "type": "deletion",
            })
        elif ":" in item:
            gene, change = item.split(":", 1)
            deletions.append({
                "gene": gene, "ref_aa": "", "codon": None, "alt_aa": "-",
                "raw": item, "hgvs_p": f"p.{change}", "type": "deletion",
            })

    insertions = []
    for item in split_csv_list(row.get("aaInsertions", "")):
        m = AA_INS_RE.match(item)
        if m:
            insertions.append({
                "gene": m.group(1),
                "codon": int(m.group(2)),
                "inserted_aa": m.group(3),
                "raw": item,
                "hgvs_p": f"p.{m.group(2)}ins{m.group(3)}",
                "type": "insertion",
            })
        elif ":" in item:
            parts = item.split(":")
            insertions.append({
                "gene": parts[0], "codon": None, "inserted_aa": ":".join(parts[1:]),
                "raw": item, "hgvs_p": f"p.{':'.join(parts[1:])}", "type": "insertion",
            })

    return substitutions, deletions, insertions


def _to_float(val, default=None):
    try:
        return float(val) if val not in (None, "") else default
    except ValueError:
        return default


def _to_int(val, default=0):
    try:
        return int(float(val)) if val not in (None, "") else default
    except ValueError:
        return default


def build_genotype_record(sample_id, consensus_path, serotype_data, coverage_data,
                          nextclade_csv_path=None, min_completeness=0.80):
    _hdr, seq = read_first_fasta(consensus_path) if consensus_path and os.path.exists(consensus_path) else ("", "")
    seq_stats = sequence_stats(seq)
    serotype = serotype_data.get("serotype", "unknown")

    nc_row = None
    if nextclade_csv_path and os.path.exists(nextclade_csv_path) and os.path.getsize(nextclade_csv_path) > 0:
        with open(nextclade_csv_path, newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter=";")
            nc_row = next(reader, None)

    qc_reasons = []
    if serotype == "unknown":
        qc_reasons.append("serotype_unknown")
    if seq_stats["length"] == 0:
        qc_reasons.append("empty_consensus")
    elif seq_stats["completeness"] < min_completeness:
        qc_reasons.append(
            f"consensus_completeness_{seq_stats['completeness']:.2%}_below_{min_completeness:.0%}"
        )

    raw_clade = "unknown"
    raw_genotype = "Unclassified"
    raw_major = None
    raw_minor = None
    nc_qc_status = "not_run"
    nc_qc_score = None
    nc_coverage = None
    aa_subs, aa_dels, aa_ins = [], [], []
    nuc_subs, nuc_dels, nuc_ins, frameshifts = [], [], [], []
    nc_errors = []
    nc_warnings = []
    total_missing = seq_stats["n_bases"]
    total_non_acgtn = seq_stats["non_acgtn_bases"]

    if nc_row is not None:
        raw_clade = (nc_row.get("clade") or "unknown").strip()
        raw_genotype, raw_major, raw_minor = parse_clade_label(raw_clade)
        nc_qc_status = (nc_row.get("qc.overallStatus") or "unknown").strip().lower()
        nc_qc_score = _to_float(nc_row.get("qc.overallScore"))
        nc_coverage = _to_float(nc_row.get("coverage"))
        aa_subs, aa_dels, aa_ins = parse_aa_mutations(nc_row)
        nuc_subs = split_csv_list(nc_row.get("substitutions", ""))
        nuc_dels = split_csv_list(nc_row.get("deletions", ""))
        nuc_ins = split_csv_list(nc_row.get("insertions", ""))
        frameshifts = split_csv_list(nc_row.get("frameShifts", ""))
        nc_errors = split_csv_list(nc_row.get("errors", ""))
        nc_warnings = split_csv_list(nc_row.get("warnings", ""))
        total_missing = _to_int(nc_row.get("totalMissing"), total_missing)
        total_non_acgtn = _to_int(nc_row.get("totalNonACGTNs"), total_non_acgtn)

        if nc_qc_status == "bad":
            qc_reasons.append(f"nextclade_qc_bad (score={nc_qc_score})")
        if nc_errors:
            qc_reasons.append(f"nextclade_errors: {'; '.join(nc_errors)}")
    else:
        if serotype != "unknown":
            qc_reasons.append("nextclade_no_result")

    qc_status = "PASS" if not qc_reasons else "FAIL"

    if qc_status == "PASS":
        genotype = raw_genotype
        major_lineage = raw_major or "Unassigned"
        minor_lineage = raw_minor or "Unassigned"
        clade = raw_clade
        if nc_qc_status == "good" and seq_stats["completeness"] >= 0.90:
            confidence = "high"
        elif nc_qc_status in ("good", "mediocre") and seq_stats["completeness"] >= min_completeness:
            confidence = "medium"
        else:
            confidence = "low"
    else:
        genotype = "Indeterminate"
        major_lineage = "Indeterminate"
        minor_lineage = "Indeterminate"
        clade = "Indeterminate"
        confidence = "none" if serotype == "unknown" or seq_stats["length"] == 0 else "low"

    return {
        "sample_id": sample_id,
        "serotype": serotype,
        "reference_accession": serotype_data.get("selected_accession") or coverage_data.get("reference"),
        "genotype": genotype,
        "major_lineage": major_lineage,
        "minor_lineage": minor_lineage,
        "clade": clade,
        "raw_clade": raw_clade,
        "raw_genotype": raw_genotype,
        "raw_major_lineage": raw_major,
        "raw_minor_lineage": raw_minor,
        "confidence": confidence,
        "method": "nextclade",
        "qc_status": qc_status,
        "qc_reasons": qc_reasons,
        "min_completeness_threshold": min_completeness,
        "sequence_stats": seq_stats,
        "nextclade": {
            "qc_status": nc_qc_status,
            "qc_score": nc_qc_score,
            "coverage": nc_coverage,
            "total_missing": total_missing,
            "total_non_acgtn": total_non_acgtn,
            "frameshifts": frameshifts,
            "warnings": nc_warnings,
            "errors": nc_errors,
            "total_nuc_substitutions": len(nuc_subs),
            "total_nuc_deletions": len(nuc_dels),
            "total_nuc_insertions": len(nuc_ins),
        },
        "aa_substitutions": aa_subs,
        "aa_deletions": aa_dels,
        "aa_insertions": aa_ins,
        "mutations": [m["raw"] for m in aa_subs],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--serotype_json", required=True)
    parser.add_argument("--coverage_json", required=True)
    parser.add_argument("--nextclade_csv", default=None)
    parser.add_argument("--min_completeness", type=float, default=0.80)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    with open(args.serotype_json, encoding="utf-8") as handle:
        serotype_data = json.load(handle)
    with open(args.coverage_json, encoding="utf-8") as handle:
        coverage_data = json.load(handle)

    record = build_genotype_record(
        sample_id=args.sample_id,
        consensus_path=args.consensus,
        serotype_data=serotype_data,
        coverage_data=coverage_data,
        nextclade_csv_path=args.nextclade_csv,
        min_completeness=args.min_completeness,
    )

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2)

    print(
        f"[{args.sample_id}] serotype={record['serotype']} clade={record['clade']} "
        f"qc={record['qc_status']} completeness={record['sequence_stats']['completeness']:.2%} "
        f"aa_sub={len(record['aa_substitutions'])} aa_del={len(record['aa_deletions'])} aa_ins={len(record['aa_insertions'])}"
    )


if __name__ == "__main__":
    main()
