#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRIMMOMATIC; CHOPPER }   from './workflows/trimming.nf'
include { HOST_REMOVAL }           from './workflows/host_removal.nf'
include { SEROTYPING }             from './workflows/serotyping.nf'
include { ILLUMINA }               from './workflows/illumina.nf'
include { NANOPORE }               from './workflows/nanopore.nf'
include { GENOTYPING }             from './workflows/genotyping.nf'
include { FHIR }                   from './workflows/fhir.nf'
include { MERGE_CLINICAL_DATA }    from './workflows/merge_clinical_data.nf'
include { VALIDATE }               from './workflows/validate_fhir.nf'
include { UPLOAD_FHIR }           from './workflows/upload_fhir.nf'
include { GENERATE_REPORT; GENERATE_DENGUE_REPORTS }        from './workflows/report.nf'
include { VERSIONS }               from './workflows/utils.nf'

workflow {

    log.info """
    Dengue Virus Analysis Pipeline (v${params.version})
    Developed by SPHERES OUCRU-ID
    Documentation: https://denv-pipeline-docs.readthedocs.io/
"""

    illumina_reads_ch = Channel
        .fromFilePairs("${params.reads_dir}/*_{1,2}_illumina.fastq.gz", checkIfExists: false)
        .map { id, files -> tuple(id, files) }

    nanopore_reads_ch = Channel
        .fromPath("${params.reads_dir}/*_ont.fastq.gz", checkIfExists: false)
        .map { file -> tuple(file.baseName.replaceFirst(/_ont$/, ''), file) }

    illumina_trimmed = TRIMMOMATIC(illumina_reads_ch)
    nanopore_trimmed = CHOPPER(nanopore_reads_ch)
    all_reads = illumina_trimmed.trimmed_reads.mix(nanopore_trimmed.trimmed_reads)
    host_removal_out = HOST_REMOVAL(all_reads, 'mixed')
    illumina_clean = host_removal_out.host_removed
        .filter { id, reads -> reads instanceof List }
    nanopore_clean = host_removal_out.host_removed
        .filter { id, reads -> !(reads instanceof List) }
    
    all_clean = illumina_clean.mix(nanopore_clean)
    serotyping_out = SEROTYPING(all_clean)
    illumina_with_serotype = illumina_clean.join(serotyping_out.serotype_info)
    nanopore_with_serotype = nanopore_clean.join(serotyping_out.serotype_info)
    illumina_out = ILLUMINA(illumina_with_serotype)
    nanopore_out = NANOPORE(nanopore_with_serotype)
    
    all_consensus = illumina_out.consensus.mix(nanopore_out.consensus)
    all_qc = illumina_out.qc_report.mix(nanopore_out.qc_report)
    all_vcfs = illumina_out.vcf.mix(nanopore_out.vcf)
    all_coverage = illumina_out.coverage.mix(nanopore_out.coverage)

    genotype_out = GENOTYPING(
        all_consensus, 
        serotyping_out.serotype_info
    )

    lineage_files = serotyping_out.serotype_info
        .join(genotype_out.genotype_info)
        .join(genotype_out.raw_nextclade)
        .map { sample_id, serotype_json, genotype_json, nextclade_json ->
            [serotype_json, genotype_json, nextclade_json]
        }
        .collect()

    fhir_input_ch = all_consensus
        .join(all_coverage)
        .map { id, cons, cov -> tuple(id, cons, cov) }

    patient_metadata_ch      = Channel.fromPath(params.patient_metadata,      checkIfExists: false).first()
    org_metadata_ch          = Channel.fromPath(params.organization_metadata, checkIfExists: false).first()
    practitioner_metadata_ch = Channel.fromPath(params.practitioner_metadata, checkIfExists: false).first()
    fhir_out = FHIR(fhir_input_ch, lineage_files, org_metadata_ch)
    merged_fhir = MERGE_CLINICAL_DATA(
        fhir_out.fhir_output,
        patient_metadata_ch,
        org_metadata_ch,
        practitioner_metadata_ch
    )
    
    validated_fhir = VALIDATE(merged_fhir.merged_fhir)

    upload_out = UPLOAD_FHIR(validated_fhir.validated_fhir)
    
    qc_files_only = all_qc.map { sample_id, files -> files }.flatten()
    GENERATE_REPORT(qc_files_only.collect())
    
    combined_data = serotyping_out.serotype_info
        .join(genotype_out.genotype_info, remainder: true)
        .join(all_consensus)
        .join(all_coverage)
        .map { sample_id, serotype_json, genotype_json, consensus, coverage ->
            tuple(sample_id, serotype_json, genotype_json ?: file("empty.json"), consensus, coverage)
        }
    
    GENERATE_DENGUE_REPORTS(combined_data)
    VERSIONS()
}
