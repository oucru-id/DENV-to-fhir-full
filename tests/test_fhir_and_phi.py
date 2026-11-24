import json
import subprocess
import sys
from pathlib import Path

from clinical_to_fhir import main as clinical_main
from denvfhir_genomic import build_genomic_bundle
from denvfhir_profiles import get_profile
from link_fhir_bundle import main as link_main


def _sample_inputs(tmp_path, sample_id="LB_001"):
    cons = tmp_path / f"{sample_id}.consensus.fasta"
    seq = "ACGT" * 250
    cons.write_text(f">{sample_id}\n{seq}\n", encoding="utf-8")
    sero = {
        "sample_id": sample_id,
        "serotype": "DENV-1",
        "confidence": "high",
        "sylvatic": False,
        "selected_accession": "NC_001477.1",
        "mapped_denv_reads": 5000,
        "primary_fraction": 0.995,
        "coinfection_suspected": False,
    }
    geno = {
        "sample_id": sample_id,
        "serotype": "DENV-1",
        "reference_accession": "NC_001477.1",
        "genotype": "Genotype V",
        "major_lineage": "Lineage D",
        "minor_lineage": "Sub-lineage 1.1",
        "clade": "1V_D.1.1",
        "confidence": "high",
        "qc_status": "PASS",
        "qc_reasons": [],
        "sequence_stats": {
            "length": len(seq), "non_n_bases": len(seq), "n_bases": 0,
            "non_acgtn_bases": 0, "completeness": 1.0, "n_percentage": 0.0, "gc_percentage": 50.0,
        },
        "aa_substitutions": [
            {"gene": "E", "ref_aa": "D", "codon": 37, "alt_aa": "N", "raw": "E:D37N",
             "hgvs_p": "p.D37N", "type": "substitution"},
        ],
        "aa_deletions": [
            {"gene": "NS1", "ref_aa": "K", "codon": 12, "alt_aa": "-", "raw": "NS1:K12-",
             "hgvs_p": "p.K12del", "type": "deletion"},
        ],
        "aa_insertions": [],
    }
    cov = {
        "sample_id": sample_id,
        "assessed": True,
        "status": "assessed",
        "reference": "NC_001477.1",
        "reference_length": len(seq),
        "min_depth_threshold": 10,
        "mean_depth": 450.2,
        "median_depth": 420.0,
        "breadth_1x": 0.995,
        "breadth_min_depth": 0.99,
        "breadth_30x": 0.98,
        "masked_bases": 10,
    }
    return cons, seq, sero, geno, cov


def _default_config():
    return {
        "pipeline_version": "2.0.0",
        "source_repository": "https://github.com/oucru-id/DENV-to-fhir-full",
        "consensus_min_depth": 10,
        "consensus_min_af": 0.5,
        "variant_min_quality": 20.0,
        "min_completeness": 0.80,
        "primer_trimmed": False,
    }


def test_genomic_bundle_generic_embeds_sequence_and_fixes_v1_bugs(tmp_path):
    cons, seq, sero, geno, cov = _sample_inputs(tmp_path, sample_id="LB_001")
    bundle = build_genomic_bundle(
        "LB_001", str(cons), sero, geno, cov,
        profile=get_profile("generic"),
        config=_default_config(),
        include_sequence=True,
        performer_org_id="99000001",
    )
    raw = json.dumps(bundle)
    # No sample ID truncation (LB_001 -> LB-001, never LB)
    assert "LB-001-genomic-report" in raw
    # No serology antibody LOINC 31343-7 and no human HGNC gene system
    assert "31343-7" not in raw
    assert "genenames.org" not in raw
    assert "diagnostic-implication" not in raw
    # PHI-free
    assert "Patient/LB-001-patient" in raw
    assert "9900000000000001" not in raw

    resources = [e["resource"] for e in bundle["entry"]]
    rtypes = [r["resourceType"] for r in resources]
    assert "DiagnosticReport" in rtypes
    assert "Device" in rtypes
    assert "Provenance" in rtypes

    cons_obs = [r for r in resources if r["resourceType"] == "Observation" and r["id"].endswith("-consensus-sequence")]
    assert len(cons_obs) == 1
    assert cons_obs[0]["valueString"] == seq

    var_obs = [r for r in resources if r["resourceType"] == "Observation" and "-var-" in r["id"]]
    assert len(var_obs) == 2  # 1 substitution + 1 deletion


def test_three_layer_link_generic_and_satusehat(tmp_path, metadata_dir):
    cons, _seq, sero, geno, cov = _sample_inputs(tmp_path, sample_id="DENV_001")
    for profile_name in ("generic", "satusehat"):
        pdir = tmp_path / profile_name
        pdir.mkdir()
        bundle = build_genomic_bundle(
            "DENV_001", str(cons), sero, geno, cov,
            profile=get_profile(profile_name),
            config=_default_config(),
            include_sequence=True,
        )
        genomic_path = pdir / "DENV_001.genomic.fhir.json"
        genomic_path.write_text(json.dumps(bundle), encoding="utf-8")

        sids = pdir / "sample_ids.txt"
        sids.write_text("DENV_001\n", encoding="utf-8")
        clin_dir = pdir / "clinical"
        clinical_main([
            "--patient_metadata", str(metadata_dir / "patient_clinical_metadata.csv"),
            "--organization_metadata", str(metadata_dir / "organization_metadata.csv"),
            "--practitioner_metadata", str(metadata_dir / "practitioner_metadata.csv"),
            "--sample_ids", str(sids),
            "--fhir_profile", profile_name,
            "--outdir", str(clin_dir),
        ])

        out_dir = pdir / "bundle"
        link_main([
            "--genomic", str(genomic_path),
            "--clinical_dir", str(clin_dir),
            "--fhir_profile", profile_name,
            "--bundle_request_method", "PUT",
            "--require_clinical",
            "--outdir", str(out_dir),
        ])

        linked_path = out_dir / "DENV-001.bundle.fhir.json"
        assert linked_path.exists()
        linked = json.loads(linked_path.read_text(encoding="utf-8"))
        assert linked["type"] == "transaction"
        rtypes = {e["resource"]["resourceType"] for e in linked["entry"]}
        assert {"Patient", "Specimen", "Encounter", "ServiceRequest", "DiagnosticReport",
                "Observation", "Organization", "Practitioner", "Device", "Provenance"} <= rtypes
        if profile_name == "satusehat":
            assert "https://fhir.kemkes.go.id/r4/StructureDefinition/Bundle" in linked["meta"]["profile"]
            assert "https://fhir.kemkes.go.id/id/nik" in json.dumps(linked)
        else:
            assert "https://fhir.kemkes.go.id/id/nik" not in json.dumps(linked)


def test_phi_guard_catches_real_nik(tmp_path, repo_root):
    bad = tmp_path / "bad.csv"
    bad.write_text("sample_id,nik\nS1," + "3171" + "010101900001\n", encoding="utf-8")
    proc = subprocess.run(
        [sys.executable, str(repo_root / "tools" / "dev" / "check_no_phi.py"), str(bad)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 1
