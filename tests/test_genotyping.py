from parse_nextclade import build_genotype_record, parse_clade_label


def _write_nc_csv(path, clade="1V_D.1.1", qc_status="good", qc_score="0",
                  aa_sub="C:V26G,E:D37N", aa_del="NS1:K12-", aa_ins="E:100:GS",
                  frameshifts=""):
    header = (
        "index;seqName;clade;qc.overallScore;qc.overallStatus;coverage;"
        "aaSubstitutions;aaDeletions;aaInsertions;substitutions;deletions;insertions;"
        "frameShifts;totalMissing;totalNonACGTNs;warnings;errors\n"
    )
    row = (
        f"0;S1;{clade};{qc_score};{qc_status};0.99;"
        f"{aa_sub};{aa_del};{aa_ins};A100G;500-502;;{frameshifts};10;0;;\n"
    )
    path.write_text(header + row, encoding="utf-8")


def test_parse_clade_labels():
    assert parse_clade_label("1V_D.1.1") == ("Genotype V", "Lineage D", "Sub-lineage 1.1")
    assert parse_clade_label("2III_C.2") == ("Genotype III", "Lineage C", "Sub-lineage 2")
    assert parse_clade_label("3III_B.3.2") == ("Genotype III", "Lineage B", "Sub-lineage 3.2")
    assert parse_clade_label("4II_B.2") == ("Genotype II", "Lineage B", "Sub-lineage 2")
    assert parse_clade_label("1I") == ("Genotype I", None, None)
    assert parse_clade_label("unknown") == ("Unclassified", None, None)


def test_genotype_pass_parses_subs_dels_ins(tmp_path):
    cons = tmp_path / "cons.fasta"
    cons.write_text(">S1\n" + ("ACGT" * 250) + "\n")
    nc = tmp_path / "nc.csv"
    _write_nc_csv(nc)
    rec = build_genotype_record(
        "S1", str(cons),
        {"serotype": "DENV-1", "selected_accession": "NC_001477.1"},
        {"assessed": True, "reference": "NC_001477.1"},
        nextclade_csv_path=str(nc),
        min_completeness=0.80,
    )
    assert rec["qc_status"] == "PASS"
    assert rec["genotype"] == "Genotype V"
    assert rec["major_lineage"] == "Lineage D"
    assert rec["minor_lineage"] == "Sub-lineage 1.1"
    assert rec["clade"] == "1V_D.1.1"
    assert rec["confidence"] == "high"
    assert len(rec["aa_substitutions"]) == 2
    assert len(rec["aa_deletions"]) == 1
    assert rec["aa_deletions"][0]["hgvs_p"] == "p.K12del"
    assert len(rec["aa_insertions"]) == 1
    assert rec["aa_insertions"][0]["hgvs_p"] == "p.100insGS"


def test_genotype_gated_indeterminate_low_completeness(tmp_path):
    cons = tmp_path / "low.fasta"
    cons.write_text(">S2\n" + ("A" * 500) + ("N" * 500) + "\n")  # 50% completeness
    nc = tmp_path / "nc.csv"
    _write_nc_csv(nc, clade="2III_C.2", qc_status="mediocre", qc_score="45")
    rec = build_genotype_record(
        "S2", str(cons),
        {"serotype": "DENV-2", "selected_accession": "NC_001474.2"},
        {"assessed": True, "reference": "NC_001474.2"},
        nextclade_csv_path=str(nc),
        min_completeness=0.80,
    )
    assert rec["qc_status"] == "FAIL"
    assert rec["genotype"] == "Indeterminate"
    assert rec["clade"] == "Indeterminate"
    assert rec["raw_clade"] == "2III_C.2"
    assert any("below_80%" in r for r in rec["qc_reasons"])


def test_genotype_gated_indeterminate_nextclade_bad(tmp_path):
    cons = tmp_path / "bad.fasta"
    cons.write_text(">S3\n" + ("ACGT" * 250) + "\n")
    nc = tmp_path / "nc.csv"
    _write_nc_csv(nc, clade="3III_B.3.2", qc_status="bad", qc_score="150")
    rec = build_genotype_record(
        "S3", str(cons),
        {"serotype": "DENV-3", "selected_accession": "NC_001475.2"},
        {"assessed": True, "reference": "NC_001475.2"},
        nextclade_csv_path=str(nc),
        min_completeness=0.80,
    )
    assert rec["qc_status"] == "FAIL"
    assert rec["genotype"] == "Indeterminate"
    assert any("nextclade_qc_bad" in r for r in rec["qc_reasons"])
