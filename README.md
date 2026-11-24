# Dengue Virus Genomic Analysis to FHIR Pipeline (DENVtoFHIR)

A Nextflow pipeline for Dengue virus (DENV-1 to DENV-4) amplicon/WGS genomic analysis from raw sequencing reads or consensus FASTAs. It performs read QC, optional host depletion, panel-based serotyping with mixed-serotype co-infection detection, amplicon primer trimming, depth-masked haploid consensus reconstruction, Nextclade genotyping/lineage classification, and exports results as HL7 FHIR R4 Genomics Reporting resources.

## Key features

- **Multi-platform input**: Illumina paired-end amplicon/WGS reads, Oxford Nanopore reads, or pre-built consensus FASTAs via `--input samplesheet.csv` (or legacy `--reads_dir`).
- **Serotyping with co-infection detection**: Maps reads against the 6-genome reference panel (`data/references/panel.tsv`, DENV-1..4 + sylvatic DENV-2/4), flags suspected mixed-serotype co-infections, and never silently falls back to DENV-1 on unmapped samples.
- **Amplicon-aware haploid consensus with low-depth masking**:
  - Optional amplicon primer trimming via `ivar trim` (`--primer_bed`).
  - Illumina haploid calling with uncapped depth (`bcftools mpileup -d 1000000 | bcftools call --ploidy 1 -mv`) and DP/QUAL/AF filters (eliminating IUPAC ambiguity inflation).
  - Nanopore variant calling via `medaka`.
  - Per-base coverage via `mosdepth` masks all positions below `--consensus_min_depth` (default `10`) with `N`, ensuring reported completeness reflects real sequencing breadth.
- **Nextclade genotyping with QC gating**: Downloads each needed serotype dataset (`community/v-gen-lab/dengue/denv1..4`) once per run (or uses pre-fetched offline datasets via `--nextclade_dir`). Samples below `--min_completeness` (default `0.80`) or with Nextclade QC `bad` are reported as `Indeterminate` with explicit reasons.
- **Three-layer FHIR R4 output** (aligned with TBtoFHIR 2.0 and KPtoFHIR 2.0):
  - `fhir/genomic/<id>.genomic.fhir.json`: PHI-free `collection` Bundle (DiagnosticReport, Classification, Region-Studied QC, embedded Consensus Sequence, Variant Observations including AA substitutions/deletions/insertions, Device, Provenance).
  - `fhir/clinical/` (with `--patient_metadata`): Patient, Specimen, Encounter, ServiceRequest, Organization, Practitioner, PractitionerRole.
  - `fhir/bundle/<id>.bundle.fhir.json` (with `--patient_metadata`): self-contained FHIR `transaction` Bundle per sample.
  - Supports `--fhir_profile generic` (default, HL7 Genomics Reporting IG 3.0.0) and `--fhir_profile satusehat` (Indonesian Ministry of Health Kemkes extensions/identifiers).
  - Optional HL7 FHIR validation (`--skip_fhir_validation false`) and gated upload (`--fhir_server_url`) using Nextflow secrets (`DENVFHIR_CLIENT_ID`, `DENVFHIR_CLIENT_SECRET`, optional `DENVFHIR_API_KEY`).

## Quick start

```bash
# Wiring check (no containers needed)
nextflow run . -profile test_stub -stub-run

# End-to-end on the 4 canonical DENV RefSeq genomes (Docker)
nextflow run . -profile test_fasta,docker

# Amplicon Illumina samplesheet (Docker)
nextflow run . -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --skip_host_removal true \
  --primer_bed primers/
```

## Samplesheet format

```csv
sample,platform,fastq_1,fastq_2,fasta
DENV001,illumina,/data/DENV001_1.fastq.gz,/data/DENV001_2.fastq.gz,
DENV002,ont,/data/DENV002.fastq.gz,,
DENV003,fasta,,,/data/DENV003.fasta
```

## Output structure

```
results/
├── pipeline_info/          execution report, timeline, trace, DAG, software_versions.yml
├── qc/
│   ├── fastqc/             FastQC reports (Illumina)
│   ├── nanoplot/           NanoStats & NanoPlot HTML (ONT)
│   ├── host_removal/       Hostile logs (when enabled)
│   ├── coverage/           <id>.coverage.json, <id>.mask.bed, mosdepth outputs
│   └── multiqc/            multiqc_report.html
├── serotyping/             <id>.serotype.json
├── variants/               <id>.filtered.vcf.gz(.tbi)
├── consensus/              <id>.consensus.fasta
├── genotyping/             <id>.nextclade.csv, <id>.nextclade.json, <id>.genotype.json
├── reports/                <id>.summary_report.txt
└── fhir/
    ├── genomic/            <id>.genomic.fhir.json           (PHI-free, includes consensus sequence)
    ├── clinical/           clinical resources, clinical_summary.tsv   (personal data)
    ├── bundle/             <id>.bundle.fhir.json, link_summary.tsv    (personal data)
    ├── validation/         validation_summary.tsv, validation_<n>.json (when validation enabled)
    └── upload/             <id>.bundle.upload.json
```
