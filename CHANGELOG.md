# Changelog

## [2.0.0]

### Architecture (aligned with TBtoFHIR 2.0 and KPtoFHIR 2.0)
- Reorganized Nextflow DSL2 codebase into `main.nf`, `workflows/denvfhir.nf`, `subworkflows/local/`, `modules/local/` (single process per module with pinned BioContainers/conda and `stub:` blocks), `conf/`, and executable scripts on `PATH` in `bin/`.
- Added `--input` samplesheet support (`illumina`, `ont`, `fasta`) validated by `assets/schema_input.json` and `nextflow_schema.json` (`nf-schema@2.5.1`). Renamed `--results_dir` to `--outdir`.
- Replaced monolithic merged FHIR generation with the 3-layer FHIR architecture: PHI-free `FHIR_GENOMIC` (with embedded consensus genome sequence via `--fhir_include_sequence`, default `true`), optional `FHIR_CLINICAL`, `FHIR_LINK`, optional batched `FHIR_VALIDATE`, and secret-authenticated `FHIR_UPLOAD` (`DENVFHIR_CLIENT_ID`, `DENVFHIR_CLIENT_SECRET`, `DENVFHIR_API_KEY`).
- Added `--fhir_profile generic` (default, HL7 Genomics Reporting IG 3.0.0) and `--fhir_profile satusehat` (opt-in Kemkes extensions/identifiers).
- Added unit tests (`pytest`), `test_stub` / `test_fasta` profiles, `tools/dev/check_no_phi.py`, `.gitleaks.toml`, and GitHub Actions CI.

### Processing & Scientific Fixes
- **Depth-masked consensus**: Added `mosdepth` per-base coverage analysis and low-depth masking (`bcftools consensus -m <mask.bed>`) below `--consensus_min_depth` (default `10`), eliminating unmasked reference bleed-through.
- **Haploid amplicon calling**: Replaced capped diploid `bcftools mpileup | bcftools call -mv` with uncapped depth (`-d 1000000`), `--ploidy 1`, `bcftools norm`, and DP/QUAL/AF filtering, plus optional `ivar trim` (`--primer_bed`) for amplicon primer trimming.
- **Nanopore variant calling**: Switched ONT consensus calling to `medaka` (`--medaka_model`) + `NanoPlot` read QC.
- **Serotyping**: Unified reference panel mapping (`data/references/panel.tsv`), fixed sylvatic DENV-2/4 reference selection, added mixed-serotype co-infection flagging, and removed silent fallback to DENV-1 on unmapped samples.
- **Nextclade genotyping**: Downloads each needed serotype dataset once per run (or stages offline via `--nextclade_dir`), captures amino acid substitutions, deletions, and insertions, and gates genotype calls as `Indeterminate` when consensus completeness `< --min_completeness` or Nextclade QC is `bad`.
- **FHIR terminology & IDs**: Removed invalid LOINC `31343-7` (antibody serology) and human HGNC codes for viral proteins; fixed sample ID truncation on underscores and lineage-file prefix collisions; added `Device` and `Provenance` resources.
