# Changelog

## [2.0.0]

- Added unit tests (`pytest`), `test_stub` / `test_fasta` profiles, `tools/dev/check_no_phi.py`, `.gitleaks.toml`, and GitHub Actions CI.
- **Depth-masked consensus**: Added `mosdepth` per-base coverage analysis and low-depth masking (`bcftools consensus -m <mask.bed>`) below `--consensus_min_depth` (default `10`).
- **Haploid amplicon calling**: Replaced capped diploid `bcftools mpileup | bcftools call -mv` with uncapped depth (`-d 1000000`), `--ploidy 1`, `bcftools norm`, and DP/QUAL/AF filtering, plus optional `ivar trim` (`--primer_bed`) for amplicon primer trimming.
- **Nanopore variant calling**: Switched ONT consensus calling to `medaka` (`--medaka_model`) + `NanoPlot` read QC.
- **Serotyping**: Unified reference panel mapping (`data/references/panel.tsv`), fixed sylvatic DENV-2/4 reference selection, added mixed-serotype co-infection flagging.
- **Nextclade genotyping**: Downloads each needed serotype dataset once per run (or stages offline via `--nextclade_dir`), captures amino acid substitutions, deletions, and insertions, and gates genotype calls as `Indeterminate` when consensus completeness `< --min_completeness` or Nextclade QC is `bad`.
- **FHIR terminology & IDs**: Removed invalid LOINC `31343-7` (antibody serology) and human HGNC codes for viral proteins; fixed sample ID truncation on underscores and lineage-file prefix collisions; added `Device` and `Provenance` resources.
