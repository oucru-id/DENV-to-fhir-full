import gzip
import json
from pathlib import Path

from serotype_classifier import classify_paf
from summarize_coverage import make_unmapped, prepare_fasta, summarize_mosdepth


def _paf_line(qname, tname, nmatch=145, alen=150, mapq=60):
    return f"{qname}\t150\t0\t150\t+\t{tname}\t10735\t100\t250\t{nmatch}\t{alen}\t{mapq}\n"


def test_serotype_clean_denv1(tmp_path, panel_tsv):
    paf = tmp_path / "s1.paf"
    paf.write_text("".join(_paf_line(f"r{i}", "DENV-1|NC_001477.1") for i in range(120)))
    res = classify_paf(str(paf), str(panel_tsv), "LB_001", "illumina", total_reads=150)
    assert res["serotype"] == "DENV-1"
    assert res["selected_reference"] == "NC_0014771.fasta"
    assert res["selected_accession"] == "NC_001477.1"
    assert res["nextclade_dataset"] == "community/v-gen-lab/dengue/denv1"
    assert res["sylvatic"] is False
    assert res["coinfection_suspected"] is False
    assert res["confidence"] == "high"


def test_serotype_sylvatic_selection(tmp_path, panel_tsv):
    paf = tmp_path / "syl.paf"
    paf.write_text("".join(_paf_line(f"r{i}", "DENV-2_SYLVATIC|EF105380.1") for i in range(50)))
    res = classify_paf(str(paf), str(panel_tsv), "SYL_02", "illumina", total_reads=60)
    assert res["serotype"] == "DENV-2"
    assert res["sylvatic"] is True
    assert res["selected_reference"] == "EF105380.fasta"
    assert res["selected_accession"] == "EF105380.1"


def test_serotype_coinfection_flagged(tmp_path, panel_tsv):
    paf = tmp_path / "mix.paf"
    lines = [_paf_line(f"a{i}", "DENV-1|NC_001477.1") for i in range(80)]
    lines += [_paf_line(f"b{i}", "DENV-3|NC_001475.2") for i in range(20)]
    paf.write_text("".join(lines))
    res = classify_paf(str(paf), str(panel_tsv), "MIX01", "illumina", total_reads=100,
                       min_reads=10, coinfection_min_fraction=0.10)
    assert res["serotype"] == "DENV-1"
    assert res["coinfection_suspected"] is True
    assert res["secondary_serotype"] == "DENV-3"
    assert res["secondary_reads"] == 20
    assert res["secondary_fraction"] == 0.2


def test_serotype_insufficient_reads_never_defaults_denv1(tmp_path, panel_tsv):
    paf = tmp_path / "low.paf"
    paf.write_text("".join(_paf_line(f"r{i}", "DENV-2|NC_001474.2") for i in range(4)))
    res = classify_paf(str(paf), str(panel_tsv), "NEG01", "illumina", total_reads=100, min_reads=10)
    assert res["serotype"] == "unknown"
    assert res["selected_reference"] == "UNKNOWN"
    assert res["nextclade_dataset"] == "unknown"
    assert res["confidence"] == "none"


def test_mosdepth_masking_and_breadth(tmp_path):
    ref = tmp_path / "ref.fasta"
    ref.write_text(">NC_001477.1\n" + ("ACGT" * 25) + "\n")  # 100 bp
    bed = tmp_path / "sample.per-base.bed.gz"
    with gzip.open(bed, "wt") as handle:
        # 0..10: depth 0 (masked), 10..20: depth 5 (masked), 20..90: depth 50 (kept), 90..100: omitted -> depth 0 (masked)
        handle.write("NC_001477.1\t0\t10\t0\n")
        handle.write("NC_001477.1\t10\t20\t5\n")
        handle.write("NC_001477.1\t20\t90\t50\n")
    mask = tmp_path / "mask.bed"
    summary = summarize_mosdepth("S1", str(bed), str(ref), min_depth=10, output_mask_bed=str(mask))

    assert summary["assessed"] is True
    assert summary["reference_length"] == 100
    assert summary["breadth_1x"] == 0.80
    assert summary["breadth_min_depth"] == 0.70
    assert summary["breadth_30x"] == 0.70
    assert summary["masked_bases"] == 30
    assert summary["mean_depth"] == 35.5
    assert summary["median_depth"] == 50.0
    mask_lines = [ln.strip() for ln in mask.read_text().splitlines() if ln.strip()]
    assert mask_lines == ["NC_001477.1\t0\t20", "NC_001477.1\t90\t100"]


def test_prepare_fasta_and_unmapped(tmp_path):
    inp = tmp_path / "in.fasta"
    inp.write_text(">orig_header desc\nACGTNNNNACGT\n")
    out = tmp_path / "out.fasta"
    res = prepare_fasta("FA_01", str(inp), str(out))
    assert res["status"] == "submitted_consensus"
    assert res["sequence_stats"]["length"] == 12
    assert res["sequence_stats"]["n_bases"] == 4
    assert res["sequence_stats"]["completeness"] == round(8 / 12, 4)
    assert out.read_text().splitlines()[0] == ">FA_01"

    empty_out = tmp_path / "unmapped.fasta"
    unm = make_unmapped("UNM_01", str(empty_out))
    assert unm["status"] == "unmapped"
    assert unm["sequence_stats"]["completeness"] == 0.0
