# Pan-Serotype Dengue Virus Genomics to FHIR Pipeline (DENVtoFHIR)

A platform-agnostic Nextflow pipeline for Dengue virus genomic analysis from raw sequencing reads,
producing HL7 FHIR R4 genomics bundles (IG v3.0.0). [Full documentation](https://denv-pipeline-docs.readthedocs.io/)

## Key Features

- **Multi-platform:** Illumina paired-end short reads and Oxford Nanopore (ONT) long reads.
- **Serotyping & Genotyping:** Pan-serotype classification (DENV-1–4) with genotype and lineage resolution via Nextclade.
- **Mutation Detection:** Amino acid substitution identification from consensus sequences.
- **Quality Control:** Per-sample FastQC reports aggregated by MultiQC.
- **FHIR Compliance:** HL7 FHIR R4 bundles with Classification, Consensus Sequence, Variant, and DiagnosticReport resources.
- **Clinical Integration:** Merges genomic results with patient, organization, and practitioner metadata.

## Installation

### Setup

```bash
git clone https://github.com/oucru-id/DENV-to-fhir-full.git
cd dengue-to-fhir-full
```

### Dependencies

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Verify
nextflow -v
```

Required tools: `minimap2`, `samtools`, `bwa-mem2`, `trimmomatic`, `chopper`, `hostile`, `nextclade`, `fastqc`, `multiqc`, `python3` (with `biopython`, `pandas`).

## Directory Structure

```
denv-to-fhir-full
├── main.nf                             # Main workflow
├── nextflow.config                     # Configuration and parameters
├── workflows/
│   ├── illumina.nf                     # Illumina sub-workflow
│   ├── nanopore.nf                     # Nanopore sub-workflow
│   ├── trimming.nf                     # Read trimming
│   ├── host_removal.nf                 # Host read removal
│   ├── serotyping.nf                   # Serotype classification
│   ├── genotyping.nf                   # Genotype classification
│   ├── fhir.nf                         # FHIR Bundle generation
│   ├── validate_fhir.nf                # FHIR validation
│   ├── merge_clinical_data.nf          # Clinical metadata merge
│   ├── upload_fhir.nf                  # FHIR server upload
│   ├── report.nf                       # QC and sample report generation
│   └── utils.nf                        # Utility functions
├── scripts/
│   ├── annotated_to_fhir.py            # Consensus-to-FHIR converter
│   ├── clinical_metadata_parser.py     # Patient/org/practitioner parser
│   ├── serotype_classification.py      # Serotype classifier
│   ├── generate_dengue_report.py       # Per-sample report generator
│   ├── get_access_token.py             # FHIR access token retriever
│   ├── upload_fhir.py                  # FHIR upload script
│   ├── merge_clinical_fhir.py          # FHIR genomics + clinical data merger
│   └── get_versions.py                 # Software version collector
├── data/
│   ├── NGS/                            # Input FASTQ files
│   ├── references/                     # Reference genomes (DENV-1–4 & Sylvatic)
│   ├── patient_clinical_metadata.csv   # Patient metadata
│   ├── organization_metadata.csv       # Organization metadata
│   └── practitioner_metadata.csv       # Practitioner metadata
└── tools/
    └── fhir-validator.jar              # HL7 FHIR validator
```

## Input Data

### Illumina Reads

Place paired-end FASTQ files in `data/NGS/`:

```
data/NGS/SAMPLE_1_illumina.fastq.gz
data/NGS/SAMPLE_2_illumina.fastq.gz
```

### Nanopore Reads

Place single-end FASTQ files in `data/NGS/`:

```
data/NGS/SAMPLE_ont.fastq.gz
```

### Metadata

Fill in the three CSV files under `data/` with appropriate metadata for patients, organizations, and practitioners.

## Usage

### Get Access Token (FHIR Upload)

```bash
python scripts/get_access_token.py
```

### Basic Run

```bash
nextflow run main.nf
```

### Run with FHIR Upload

> Get the access token first before running with upload.

```bash
nextflow run main.nf \
  --fhir_server_url "https://<BASE_URL>/fhir"
```

## Output Structure

```
results/
├── qc/
│   └── multiqc_report.html             # Aggregated QC report
├── consensus/
│   └── *.consensus.fasta               # Per-sample consensus sequences
├── serotyping/
│   └── *.serotype.json                 # Per-sample serotype results
├── genotyping/
│   └── *.genotype_lineage.json         # Per-sample genotype/lineage results
├── fhir/
│   └── *.fhir.json                     # FHIR genomics bundles
├── fhir_merged/
│   └── *.merged.fhir.json              # FHIR bundles with clinical data
├── fhir_validated/
│   └── *.validation.txt                # FHIR validation results
├── fhir_upload/
│   └── *.upload.json                   # FHIR upload results
├── reports/
│   └── *.summary_report.txt            # Per-sample summary reports
├── runningstat/
│   ├── execution.html                  # Nextflow execution report
│   ├── timeline.html                   # Timeline report
│   └── dag.html                        # Workflow DAG
└── software_versions.yml               # Software version manifest
```

## Support

[GitHub Issues](https://github.com/oucru-id/DENV-to-fhir-full/issues)
