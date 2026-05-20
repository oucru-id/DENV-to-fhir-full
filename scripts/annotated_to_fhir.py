#!/usr/bin/env python3

import argparse
import json
import os
import sys
import uuid
import gzip
from datetime import datetime, timezone
from Bio import SeqIO

sys.path.insert(0, os.path.dirname(__file__))
from clinical_metadata_parser import load_organization_metadata

def main():
    parser = argparse.ArgumentParser(description='Convert Dengue analysis results to FHIR observations')
    parser.add_argument('--input', required=True, help='Consensus FASTA file')
    parser.add_argument('--output', required=True, help='Output FHIR JSON file')
    parser.add_argument('--lineage_dir', required=True, help='Directory with lineage JSON files')
    parser.add_argument('--coverage_file', required=False, help='samtools coverage output file')
    parser.add_argument('--organization_metadata', required=False, help='Path to organization metadata CSV file')
    
    args = parser.parse_args()

    org_data = {}
    if args.organization_metadata and os.path.exists(args.organization_metadata):
        org_data = load_organization_metadata(args.organization_metadata)
    
    sample_id = extract_sample_id(args.input)
    seq_data = read_consensus(args.input)
    coverage_stats = read_coverage_stats(args.coverage_file)
    
    serotype_data, genotype_data, nextclade_data = read_lineage_data(args.lineage_dir, sample_id)
    fhir_bundle = create_fhir_observations(sample_id, seq_data, serotype_data, genotype_data, nextclade_data, coverage_stats, org_data)
    
    with open(args.output, 'w') as f:
        json.dump(fhir_bundle, f, indent=2)

def extract_sample_id(filename):
    basename = os.path.basename(filename)
    return basename.split('.')[0].split('_')[0]

def read_consensus(fasta_file):
    try:
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                sequence = str(record.seq)
                n_count = sequence.count('N')
                return {
                    'sequence': sequence,
                    'length': len(sequence),
                    'n_percentage': (n_count / len(sequence)) * 100 if len(sequence) > 0 else 0,
                    'gc_percentage': ((sequence.count('G') + sequence.count('C')) / len(sequence)) * 100 if len(sequence) > 0 else 0
                }
    except:
        return {'sequence': '', 'length': 0, 'n_percentage': 0, 'gc_percentage': 0}

def read_coverage_stats(coverage_file):
    try:
        if not coverage_file or not os.path.exists(coverage_file):
            return {"breadth": None, "mean_depth": None}
        with open(coverage_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    return {
                        "breadth":    float(parts[5]),
                        "mean_depth": float(parts[6])
                    }
        return {"breadth": None, "mean_depth": None}
    except Exception:
        return {"breadth": None, "mean_depth": None}


def read_lineage_data(lineage_dir, sample_id):
    serotype_data = {}
    genotype_data = {}
    nextclade_mutations = []
    
    for filename in os.listdir(lineage_dir):
        if sample_id in filename:
            filepath = os.path.join(lineage_dir, filename)
            try:
                if filename.endswith('.serotype.json'):
                    with open(filepath, 'r') as f: serotype_data = json.load(f)
                elif filename.endswith('.genotype_lineage.json'):
                    with open(filepath, 'r') as f: genotype_data = json.load(f)
                elif filename.endswith('.nextclade.csv'):
                    import csv
                    with open(filepath, 'r') as f:
                        reader = csv.DictReader(f, delimiter=';')
                        for row in reader:
                            aa_subs = row.get('aaSubstitutions', '')
                            if aa_subs:
                                for item in aa_subs.split(','):
                                    if ':' in item:
                                        gene, change = item.split(':')
                                        if len(change) >= 2:
                                            ref = change[0]
                                            alt = change[-1]
                                            pos = change[1:-1]
                                            nextclade_mutations.append({
                                                'gene': gene,
                                                'refAA': ref,
                                                'codon': pos,
                                                'queryAA': alt
                                            })
            except Exception as e:
                print(f"Error reading {filename}: {e}")
                continue
                
    return serotype_data, genotype_data, nextclade_mutations

def create_fhir_observations(sample_id, seq_data, serotype_data, genotype_data, nextclade_mutations, coverage_stats=None, org_data=None):
    org_data = org_data or {}
    org_id = org_data.get('org_id', 'unknown-org')
    timestamp = datetime.now(timezone.utc).isoformat()
    bundle = {
        "resourceType": "Bundle",
        "type": "collection",
        "timestamp": timestamp,
        "entry": []
    }
    
    serotype = serotype_data.get('serotype', 'Unknown')
    genotype = genotype_data.get('genotype', 'Unknown')
    major_lineage = genotype_data.get('major_lineage', 'Unknown')
    minor_lineage = genotype_data.get('minor_lineage', 'Unknown')
    confidence = genotype_data.get('confidence', 'low') 
    
    classification_obs = {
        "resourceType": "Observation",
        "id": f"{sample_id}-classification",
        "status": "final",
        "code": {
            "coding": [{"system": "http://loinc.org", "code": "31343-7", "display": "Dengue virus Ab [Presence] in Specimen"}],
            "text": "Dengue Virus Classification"
        },
        "valueCodeableConcept": {
            "text": f"Dengue Virus {serotype} {genotype}"
        },
        "component": [
            {
                "code": {"text": "Serotype"},
                "valueCodeableConcept": {"text": serotype}
            },
            {
                "code": {"text": "Genotype"},
                "valueCodeableConcept": {"text": genotype}
            },
            {
                "code": {"text": "Major Lineage"},
                "valueCodeableConcept": {"text": major_lineage}
            },
            {
                "code": {"text": "Minor Lineage"},
                "valueCodeableConcept": {"text": minor_lineage}
            },
            {
                "code": {"text": "Confidence"},
                "valueCodeableConcept": {"text": confidence}
            }
        ],
        "subject": {"reference": f"Patient/{sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{sample_id}-specimen"},
        "effectiveDateTime": timestamp,
        "performer": [{"reference": f"Organization/{org_id}"}]
    }
    bundle["entry"].append({"fullUrl": f"urn:uuid:{uuid.uuid4()}", "resource": classification_obs})

    if seq_data.get('sequence'):
        consensus_obs = {
            "resourceType": "Observation",
            "id": f"{sample_id}-consensus",
            "meta": {
                "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"]
            },
            "status": "final",
            "category": [
                {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory",
                            "display": "Laboratory"
                        }
                    ]
                }
            ],
            "code": {
                "coding": [
                    {
                        "system": "http://loinc.org",
                        "code": "86206-0", 
                        "display": "Whole genome sequence analysis in Blood or Tissue by Molecular genetics method"
                    }
                ],
                "text": "Viral Consensus Genome Sequence"
            },
            "valueString": seq_data['sequence'],
            "component": [
                {
                    "code": {"text": "Sequence Length"},
                    "valueQuantity": {
                        "value": seq_data['length'],
                        "unit": "bp",
                        "system": "http://unitsofmeasure.org",
                        "code": "bp"
                    }
                },
                {
                    "code": {"text": "Sequence Completeness"},
                    "valueQuantity": {
                        "value": round(100.0 - seq_data['n_percentage'], 2),
                        "unit": "%",
                        "system": "http://unitsofmeasure.org",
                        "code": "%"
                    }
                },
                {
                    "code": {"text": "Ambiguous bases (N content)"},
                    "valueQuantity": {
                        "value": round(seq_data['n_percentage'], 2),
                        "unit": "%",
                        "system": "http://unitsofmeasure.org",
                        "code": "%"
                    }
                },
                {
                    "code": {"text": "GC Content"},
                    "valueQuantity": {
                        "value": round(seq_data['gc_percentage'], 2),
                        "unit": "%",
                        "system": "http://unitsofmeasure.org",
                        "code": "%"
                    }
                }
            ] + ([
                {
                    "code": {"text": "Genome Coverage"},
                    "valueQuantity": {
                        "value": round(coverage_stats['breadth'], 2),
                        "unit": "%",
                        "system": "http://unitsofmeasure.org",
                        "code": "%"
                    }
                },
                {
                    "code": {"text": "Mean Sequencing Depth"},
                    "valueQuantity": {
                        "value": round(coverage_stats['mean_depth'], 2),
                        "unit": "x",
                        "system": "http://unitsofmeasure.org",
                        "code": "{fold}"
                    }
                }
            ] if coverage_stats and coverage_stats.get('breadth') is not None else []),
            "subject": {"reference": f"Patient/{sample_id}-patient"},
            "specimen": {"reference": f"Specimen/{sample_id}-specimen"},
            "effectiveDateTime": timestamp,
            "performer": [{"reference": f"Organization/{org_id}"}]
        }
        bundle["entry"].append({"fullUrl": f"urn:uuid:{uuid.uuid4()}", "resource": consensus_obs})

    for i, sub in enumerate(nextclade_mutations, 1):
        gene = sub.get('gene', 'Unknown')
        ref = sub.get('refAA', '')
        pos = sub.get('codon', '')
        alt = sub.get('queryAA', '')
        
        hgvs_notation = f"p.{ref}{pos}{alt}"
        
        variant_obs = {
            "resourceType": "Observation",
            "id": f"{sample_id}-obs-{i}", 
            "meta": {
                "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"]
            },
            "status": "final",
            "category": [
                {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory",
                            "display": "Laboratory"
                        }
                    ]
                },
                {
                    "coding": [
                        {
                            "system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                            "code": "diagnostic-implication",
                            "display": "Diagnostic Implication"
                        }
                    ]
                }
            ],
            "code": {
                "coding": [{"system": "http://loinc.org", "code": "69548-6", "display": "Genetic variant assessment"}]
            },
            "valueCodeableConcept": {
                "coding": [{"system": "http://loinc.org", "code": "LA9633-4", "display": "Present"}]
            },
            "component": [
                {
                    "code": {
                        "coding": [{"system": "http://loinc.org", "code": "48018-6", "display": "Gene studied [ID]"}]
                    },
                    "valueCodeableConcept": {
                        "coding": [
                            {
                                "system": "http://www.genenames.org/geneId",
                                "code": gene,
                                "display": gene
                            }
                        ],
                        "text": gene
                    }
                },
                {
                    "code": {
                        "coding": [{"system": "http://loinc.org", "code": "48005-3", "display": "Amino acid change (pHGVS)"}]
                    },
                    "valueCodeableConcept": {
                        "coding": [
                          {
                            "system": "https://varnomen.hgvs.org",
                            "code": hgvs_notation,
                            "display": hgvs_notation
                          }
                        ]
                    }
                }
            ],
            "subject": {"reference": f"Patient/{sample_id}-patient"},
            "specimen": {"reference": f"Specimen/{sample_id}-specimen"},
            "effectiveDateTime": timestamp,
            "performer": [{"reference": f"Organization/{org_id}"}]
        }
        bundle["entry"].append({"fullUrl": f"urn:uuid:{uuid.uuid4()}", "resource": variant_obs})

    return bundle

if __name__ == "__main__":
    main()
