#!/usr/bin/env python3
"""Build a PHI-free HL7 FHIR R4 genomic collection Bundle for a Dengue virus sample."""

import argparse
import base64
import json
import os

from denvfhir_common import (SOURCE_REPOSITORY, bundle_id, fhir_id, full_url, now_iso,
                             point_references_at_entries, today_compact, xhtml_text)
from denvfhir_profiles import get_profile
from denvfhir_terminology import (CS_PIPELINE_PROPERTY, HGVS, LOINC, SEQUENCE_ONTOLOGY, UCUM,
                                  VARIANT_TYPE_SO)
from summarize_coverage import read_first_fasta, sequence_stats

LAB_CATEGORY = [{
    "coding": [{
        "system": "http://terminology.hl7.org/CodeSystem/observation-category",
        "code": "laboratory",
        "display": "Laboratory",
    }]
}]


class Context:
    def __init__(self, sample_id, fid, profile, timestamp, performer_org_id=None, performer_org_name=None):
        self.sample_id = sample_id
        self.fid = fid
        self.profile = profile
        self.timestamp = timestamp
        self.performer_org_id = performer_org_id
        self.performer_org_name = performer_org_name

    @property
    def performer(self):
        if not self.performer_org_id:
            return None
        ref = {"reference": f"Organization/{self.performer_org_id}"}
        if self.performer_org_name:
            ref["display"] = self.performer_org_name
        return [ref]


def _attach_common_obs(obs, ctx, kind):
    meta = ctx.profile.meta(kind)
    if meta:
        obs["meta"] = meta
    obs["status"] = "final"
    obs["category"] = LAB_CATEGORY
    obs["subject"] = {"reference": f"Patient/{ctx.fid}-patient"}
    obs["specimen"] = {"reference": f"Specimen/{ctx.fid}-specimen"}
    obs["effectiveDateTime"] = ctx.timestamp
    if ctx.performer:
        obs["performer"] = ctx.performer
    return obs


def build_classification_observation(ctx, serotype_data, genotype_data):
    serotype = serotype_data.get("serotype", "unknown")
    genotype = genotype_data.get("genotype", "Indeterminate")
    major = genotype_data.get("major_lineage", "Indeterminate")
    minor = genotype_data.get("minor_lineage", "Indeterminate")
    clade = genotype_data.get("clade", "Indeterminate")
    qc_status = genotype_data.get("qc_status", "FAIL")
    confidence = genotype_data.get("confidence", "none")

    if serotype == "unknown":
        summary_text = "Dengue virus serotype unknown"
    elif qc_status == "PASS":
        summary_text = f"Dengue virus {serotype} ({genotype}, {clade})"
    else:
        summary_text = f"Dengue virus {serotype} (Genotype Indeterminate)"

    components = [
        {
            "code": ctx.profile.component_code("serotype"),
            "valueCodeableConcept": ctx.profile.serotype_concept(serotype),
        },
        {
            "code": ctx.profile.component_code("genotype"),
            "valueCodeableConcept": {"text": genotype},
        },
        {
            "code": ctx.profile.component_code("major-lineage"),
            "valueCodeableConcept": {"text": major},
        },
        {
            "code": ctx.profile.component_code("minor-lineage"),
            "valueCodeableConcept": {"text": minor},
        },
        {
            "code": ctx.profile.component_code("clade"),
            "valueCodeableConcept": ctx.profile.clade_concept(clade),
        },
        {
            "code": ctx.profile.component_code("classification-confidence"),
            "valueCodeableConcept": {"text": confidence},
        },
        {
            "code": ctx.profile.component_code("sequencing-qc"),
            "valueCodeableConcept": {"text": qc_status},
        },
        {
            "code": ctx.profile.component_code("sylvatic-flag"),
            "valueBoolean": bool(serotype_data.get("sylvatic", False)),
        },
        {
            "code": ctx.profile.component_code("coinfection-flag"),
            "valueBoolean": bool(serotype_data.get("coinfection_suspected", False)),
        },
    ]

    primary_frac = serotype_data.get("primary_fraction")
    if primary_frac is not None and serotype_data.get("mapped_denv_reads", 0) > 0:
        components.append({
            "code": ctx.profile.component_code("serotype-read-fraction"),
            "valueQuantity": {
                "value": round(float(primary_frac) * 100.0, 2),
                "unit": "%",
                "system": UCUM,
                "code": "%",
            },
        })

    obs = {
        "resourceType": "Observation",
        "id": fhir_id(ctx.fid, "classification"),
        "code": ctx.profile.observation_code("dengue-classification", text="Dengue Virus Classification"),
        "valueCodeableConcept": ctx.profile.serotype_concept(serotype, text=summary_text),
        "component": components,
    }
    if serotype_data.get("coinfection_suspected"):
        obs["note"] = [{
            "text": (
                f"Possible mixed-serotype co-infection: secondary serotype "
                f"{serotype_data.get('secondary_serotype')} supported by "
                f"{serotype_data.get('secondary_reads')} reads "
                f"({(serotype_data.get('secondary_fraction') or 0) * 100:.1f}%)."
            )
        }]
    return _attach_common_obs(obs, ctx, "classification")


def build_region_studied_observation(ctx, seq_stats, coverage_data, genotype_data):
    components = [
        {
            "code": ctx.profile.component_code("sequence-length"),
            "valueQuantity": {
                "value": seq_stats["length"],
                "unit": "bp",
                "system": UCUM,
                "code": "{bp}",
            },
        },
        {
            "code": ctx.profile.component_code("sequence-completeness"),
            "valueQuantity": {
                "value": round(seq_stats["completeness"] * 100.0, 2),
                "unit": "%",
                "system": UCUM,
                "code": "%",
            },
        },
        {
            "code": ctx.profile.component_code("ambiguous-bases"),
            "valueQuantity": {
                "value": round(seq_stats["n_percentage"], 2),
                "unit": "%",
                "system": UCUM,
                "code": "%",
            },
        },
        {
            "code": ctx.profile.component_code("gc-content"),
            "valueQuantity": {
                "value": round(seq_stats["gc_percentage"], 2),
                "unit": "%",
                "system": UCUM,
                "code": "%",
            },
        },
        {
            "code": ctx.profile.component_code("sequencing-qc"),
            "valueCodeableConcept": {"text": genotype_data.get("qc_status", "FAIL")},
        },
    ]

    ref_acc = coverage_data.get("reference") or genotype_data.get("reference_accession")
    if ref_acc:
        components.insert(0, {
            "code": {"coding": [{"system": LOINC, "code": "48013-7",
                                 "display": "Genomic reference sequence [ID]"}]},
            "valueCodeableConcept": {"text": ref_acc},
        })

    if coverage_data.get("assessed"):
        min_dp = coverage_data.get("min_depth_threshold", 10)
        components.extend([
            {
                "code": ctx.profile.component_code("coverage-breadth"),
                "valueQuantity": {
                    "value": round((coverage_data.get("breadth_min_depth") or 0.0) * 100.0, 2),
                    "unit": "%",
                    "system": UCUM,
                    "code": "%",
                },
            },
            {
                "code": {"coding": [{"system": LOINC, "code": "82121-5", "display": "Allelic read depth"}],
                         "text": "Mean sequencing depth"},
                "valueQuantity": {
                    "value": round(float(coverage_data.get("mean_depth") or 0.0), 2),
                    "unit": "x",
                    "system": UCUM,
                    "code": "{fold}",
                },
            },
        ])
        note_text = (
            f"Reference {ref_acc}: breadth >=1x {(coverage_data.get('breadth_1x') or 0)*100:.2f}%, "
            f"breadth >={min_dp}x {(coverage_data.get('breadth_min_depth') or 0)*100:.2f}%, "
            f"mean depth {coverage_data.get('mean_depth', 0):.1f}x, "
            f"median depth {coverage_data.get('median_depth', 0):.1f}x, "
            f"masked bases {coverage_data.get('masked_bases', 0)}."
        )
    else:
        note_text = f"Read coverage status: {coverage_data.get('status', 'not_assessed')}."

    if genotype_data.get("qc_reasons"):
        note_text += " QC flags: " + "; ".join(genotype_data["qc_reasons"])

    obs = {
        "resourceType": "Observation",
        "id": fhir_id(ctx.fid, "region-studied"),
        "code": {
            "coding": [{"system": LOINC, "code": "53041-0",
                        "display": "DNA region of interest panel"}],
            "text": "Dengue Genome Coverage and Consensus QC",
        },
        "component": components,
        "note": [{"text": note_text}],
    }
    return _attach_common_obs(obs, ctx, "region_studied")


def build_consensus_observation(ctx, seq, seq_stats):
    if not seq:
        return None
    obs = {
        "resourceType": "Observation",
        "id": fhir_id(ctx.fid, "consensus-sequence"),
        "code": ctx.profile.observation_code(
            "viral-consensus-sequence",
            extra_codings=[{
                "system": LOINC,
                "code": "86206-0",
                "display": "Whole genome sequence analysis in Blood or Tissue by Molecular genetics method",
            }],
            text="Viral Consensus Genome Sequence",
        ),
        "valueString": seq,
        "component": [
            {
                "code": ctx.profile.component_code("sequence-length"),
                "valueQuantity": {"value": seq_stats["length"], "unit": "bp", "system": UCUM, "code": "{bp}"},
            },
            {
                "code": ctx.profile.component_code("sequence-completeness"),
                "valueQuantity": {"value": round(seq_stats["completeness"] * 100.0, 2), "unit": "%",
                                  "system": UCUM, "code": "%"},
            },
            {
                "code": ctx.profile.component_code("ambiguous-bases"),
                "valueQuantity": {"value": round(seq_stats["n_percentage"], 2), "unit": "%",
                                  "system": UCUM, "code": "%"},
            },
            {
                "code": ctx.profile.component_code("gc-content"),
                "valueQuantity": {"value": round(seq_stats["gc_percentage"], 2), "unit": "%",
                                  "system": UCUM, "code": "%"},
            },
        ],
    }
    return _attach_common_obs(obs, ctx, "consensus")


def build_variant_observations(ctx, genotype_data):
    ref_acc = genotype_data.get("reference_accession")
    all_mutations = (
        list(genotype_data.get("aa_substitutions", []))
        + list(genotype_data.get("aa_deletions", []))
        + list(genotype_data.get("aa_insertions", []))
    )

    observations = []
    for idx, mut in enumerate(all_mutations, 1):
        gene = mut.get("gene", "Unknown")
        hgvs_p = mut.get("hgvs_p") or f"p.{mut.get('raw', '')}"
        mtype = mut.get("type", "substitution")
        so = VARIANT_TYPE_SO.get(mtype)

        components = [
            {
                "code": {"coding": [{"system": LOINC, "code": "48018-6", "display": "Gene studied [ID]"}]},
                "valueCodeableConcept": ctx.profile.gene_concept(gene),
            },
            {
                "code": {"coding": [{"system": LOINC, "code": "48005-3",
                                     "display": "Amino acid change (pHGVS)"}]},
                "valueCodeableConcept": {
                    "coding": [{"system": HGVS, "code": hgvs_p, "display": hgvs_p}],
                    "text": f"{ref_acc}:{gene}:{hgvs_p}" if ref_acc else f"{gene}:{hgvs_p}",
                },
            },
        ]
        if ref_acc:
            components.append({
                "code": {"coding": [{"system": LOINC, "code": "48013-7",
                                     "display": "Genomic reference sequence [ID]"}]},
                "valueCodeableConcept": {"text": ref_acc},
            })
        if so:
            components.append({
                "code": {"coding": [{"system": LOINC, "code": "48019-4", "display": "DNA change type"}]},
                "valueCodeableConcept": {
                    "coding": [{"system": SEQUENCE_ONTOLOGY, "code": so["code"], "display": so["display"]}],
                    "text": mtype,
                },
            })

        obs = {
            "resourceType": "Observation",
            "id": fhir_id(ctx.fid, "var", str(idx)),
            "code": {"coding": [{"system": LOINC, "code": "69548-6",
                                 "display": "Genetic variant assessment"}]},
            "valueCodeableConcept": {"coding": [{"system": LOINC, "code": "LA9633-4", "display": "Present"}]},
            "component": components,
        }
        observations.append(_attach_common_obs(obs, ctx, "variant"))
    return observations


def build_diagnostic_report(ctx, observations, serotype_data, genotype_data, seq_stats):
    serotype = serotype_data.get("serotype", "unknown")
    genotype = genotype_data.get("genotype", "Indeterminate")
    major = genotype_data.get("major_lineage", "Indeterminate")
    minor = genotype_data.get("minor_lineage", "Indeterminate")
    clade = genotype_data.get("clade", "Indeterminate")
    qc_status = genotype_data.get("qc_status", "FAIL")
    confidence = genotype_data.get("confidence", "none")

    if serotype == "unknown":
        conclusion = "Dengue virus serotype could not be determined (insufficient viral reads)."
    elif qc_status == "PASS":
        conclusion = (
            f"Dengue virus serotype {serotype} detected. Genotype: {genotype}. "
            f"Lineage: {major} / {minor} (clade {clade}). "
            f"Consensus completeness: {seq_stats['completeness'] * 100:.1f}%. "
            f"Confidence: {confidence.upper()}."
        )
    else:
        reasons = "; ".join(genotype_data.get("qc_reasons", []))
        conclusion = (
            f"Dengue virus serotype {serotype} detected, but genotype/lineage is Indeterminate "
            f"(QC FAIL: {reasons})."
        )
    if serotype_data.get("coinfection_suspected"):
        conclusion += (
            f" Warning: possible mixed-serotype co-infection with "
            f"{serotype_data.get('secondary_serotype')} "
            f"({(serotype_data.get('secondary_fraction') or 0) * 100:.1f}% of mapped reads)."
        )

    html = (
        '<div xmlns="http://www.w3.org/1999/xhtml">'
        "<h1>Dengue Virus Genomic Analysis Report</h1>"
        f"<p><strong>Sample ID:</strong> {xhtml_text(ctx.sample_id)}</p>"
        f"<p><strong>Report Date:</strong> {xhtml_text(ctx.timestamp)}</p>"
        f"<p><strong>Serotype:</strong> {xhtml_text(serotype)}</p>"
        f"<p><strong>Genotype:</strong> {xhtml_text(genotype)}</p>"
        f"<p><strong>Lineage / Clade:</strong> {xhtml_text(major)} / {xhtml_text(minor)} ({xhtml_text(clade)})</p>"
        f"<p><strong>QC Status:</strong> {xhtml_text(qc_status)}</p>"
        f"<p><strong>Completeness:</strong> {seq_stats['completeness'] * 100:.2f}%</p>"
        f"<p><strong>Conclusion:</strong> {xhtml_text(conclusion)}</p>"
        "</div>"
    )

    conclusion_codes = [ctx.profile.serotype_concept(serotype)]
    if qc_status == "PASS" and clade not in ("unknown", "Indeterminate"):
        conclusion_codes.append(ctx.profile.clade_concept(clade))

    report = {
        "resourceType": "DiagnosticReport",
        "id": fhir_id(ctx.fid, "genomic-report"),
        "status": "final",
        "category": [{
            "coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                        "code": "GE", "display": "Genetics"}]
        }],
        "code": ctx.profile.report_code(),
        "subject": {"reference": f"Patient/{ctx.fid}-patient"},
        "effectiveDateTime": ctx.timestamp,
        "issued": ctx.timestamp,
        "specimen": [{"reference": f"Specimen/{ctx.fid}-specimen"}],
        "result": [{"reference": f"Observation/{obs['id']}"} for obs in observations],
        "conclusion": conclusion,
        "conclusionCode": conclusion_codes,
        "presentedForm": [{
            "contentType": "text/html",
            "language": "en-US",
            "title": "Dengue Virus Genomic Analysis Report",
            "data": base64.b64encode(html.encode("utf-8")).decode("ascii"),
        }],
    }
    meta = ctx.profile.meta("report")
    if meta:
        report["meta"] = meta
    if ctx.performer:
        report["performer"] = ctx.performer
    identifier = ctx.profile.identifier(
        "diagnostic-report", f"DENV-GEN-{ctx.fid}-{today_compact()}", ctx.performer_org_id
    )
    if identifier:
        report["identifier"] = [identifier]
    return report


def build_device(profile, config):
    def prop(code, display, value):
        return {
            "type": {"coding": [{"system": CS_PIPELINE_PROPERTY, "code": code, "display": display}]},
            "valueCode": [{"text": str(value)}],
        }

    device = {
        "resourceType": "Device",
        "id": "denvfhir-pipeline",
        "status": "active",
        "deviceName": [{"name": "DENVtoFHIR", "type": "user-friendly-name"}],
        "version": [{"value": str(config["pipeline_version"])}],
        "property": [
            prop("source-repository", "Source repository", config["source_repository"]),
            prop("consensus-min-depth", "Minimum depth for consensus masking", config["consensus_min_depth"]),
            prop("consensus-min-af", "Minimum allele fraction for haploid consensus calling", config["consensus_min_af"]),
            prop("variant-min-quality", "Minimum variant quality (QUAL)", config["variant_min_quality"]),
            prop("min-completeness", "Minimum consensus completeness for genotyping", config["min_completeness"]),
            prop("primer-trimmed", "Amplicon primer trimming applied", str(config["primer_trimmed"]).lower()),
        ],
    }
    meta = profile.meta("device")
    if meta:
        device["meta"] = meta
    return device


def build_provenance(ctx, targets, device_id, description):
    agents = [{
        "type": {"coding": [{
            "system": "http://terminology.hl7.org/CodeSystem/provenance-participant-type",
            "code": "assembler",
            "display": "Assembler",
        }]},
        "who": {"reference": f"Device/{device_id}", "display": "DENVtoFHIR"},
    }]
    if ctx.performer:
        agents.append({
            "type": {"coding": [{
                "system": "http://terminology.hl7.org/CodeSystem/provenance-participant-type",
                "code": "performer",
                "display": "Performer",
            }]},
            "who": dict(ctx.performer[0]),
        })
    prov = {
        "resourceType": "Provenance",
        "id": fhir_id(ctx.fid, "provenance"),
        "target": targets,
        "recorded": ctx.timestamp,
        "activity": {
            "coding": [{"system": "http://terminology.hl7.org/CodeSystem/v3-DataOperation",
                        "code": "CREATE", "display": "create"}],
            "text": description,
        },
        "agent": agents,
    }
    meta = ctx.profile.meta("provenance")
    if meta:
        prov["meta"] = meta
    return prov


def build_genomic_bundle(sample_id, consensus_path, serotype_data, genotype_data, coverage_data,
                         profile, config, include_sequence=True,
                         performer_org_id=None, performer_org_name=None):
    fid = fhir_id(sample_id)
    timestamp = now_iso()
    ctx = Context(sample_id, fid, profile, timestamp, performer_org_id, performer_org_name)

    _hdr, seq = read_first_fasta(consensus_path) if consensus_path and os.path.exists(consensus_path) else ("", "")
    seq_stats = genotype_data.get("sequence_stats") or sequence_stats(seq)

    classification_obs = build_classification_observation(ctx, serotype_data, genotype_data)
    region_obs = build_region_studied_observation(ctx, seq_stats, coverage_data, genotype_data)
    observations = [classification_obs, region_obs]

    if include_sequence:
        cons_obs = build_consensus_observation(ctx, seq, seq_stats)
        if cons_obs:
            observations.append(cons_obs)

    variant_obs = build_variant_observations(ctx, genotype_data)
    observations.extend(variant_obs)

    report = build_diagnostic_report(ctx, observations, serotype_data, genotype_data, seq_stats)
    device = build_device(profile, config)
    targets = [{"reference": f"DiagnosticReport/{report['id']}"}] + [
        {"reference": f"Observation/{obs['id']}"} for obs in observations
    ]
    provenance = build_provenance(
        ctx, targets, device["id"],
        f"Dengue genomic analysis for {fid} ({os.path.basename(consensus_path or 'consensus')})"
    )

    resources = point_references_at_entries([report] + observations + [device, provenance])
    return {
        "resourceType": "Bundle",
        "id": bundle_id("genomic", fid),
        "type": "collection",
        "timestamp": timestamp,
        "entry": [{"fullUrl": full_url(r["resourceType"], r["id"]), "resource": r} for r in resources],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--serotype_json", required=True)
    parser.add_argument("--genotype_json", required=True)
    parser.add_argument("--coverage_json", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fhir_profile", default="generic", choices=["generic", "satusehat"])
    parser.add_argument("--include_sequence", action="store_true")
    parser.add_argument("--primer_trimmed", action="store_true")
    parser.add_argument("--performer_org_id", default=None)
    parser.add_argument("--performer_org_name", default=None)
    parser.add_argument("--identifier_system_base", default=None)
    parser.add_argument("--pipeline_version", default="2.0.0")
    parser.add_argument("--source_repository", default=SOURCE_REPOSITORY)
    parser.add_argument("--consensus_min_depth", type=int, default=10)
    parser.add_argument("--consensus_min_af", type=float, default=0.5)
    parser.add_argument("--variant_min_quality", type=float, default=20.0)
    parser.add_argument("--min_completeness", type=float, default=0.80)
    args = parser.parse_args(argv)

    with open(args.serotype_json, encoding="utf-8") as handle:
        serotype_data = json.load(handle)
    with open(args.genotype_json, encoding="utf-8") as handle:
        genotype_data = json.load(handle)
    with open(args.coverage_json, encoding="utf-8") as handle:
        coverage_data = json.load(handle)

    profile = get_profile(args.fhir_profile, args.identifier_system_base)
    config = {
        "pipeline_version": args.pipeline_version,
        "source_repository": args.source_repository,
        "consensus_min_depth": args.consensus_min_depth,
        "consensus_min_af": args.consensus_min_af,
        "variant_min_quality": args.variant_min_quality,
        "min_completeness": args.min_completeness,
        "primer_trimmed": args.primer_trimmed,
    }

    bundle = build_genomic_bundle(
        sample_id=args.sample_id,
        consensus_path=args.consensus,
        serotype_data=serotype_data,
        genotype_data=genotype_data,
        coverage_data=coverage_data,
        profile=profile,
        config=config,
        include_sequence=args.include_sequence,
        performer_org_id=args.performer_org_id,
        performer_org_name=args.performer_org_name,
    )

    with open(args.output, "w", encoding="utf-8") as out:
        json.dump(bundle, out, indent=2, ensure_ascii=False)

    print(
        f"[{args.sample_id}] genomic FHIR bundle written ({len(bundle['entry'])} entries, "
        f"profile={profile.name}, sequence_embedded={args.include_sequence})"
    )


if __name__ == "__main__":
    main()
