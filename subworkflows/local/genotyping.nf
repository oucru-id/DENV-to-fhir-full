include { NEXTCLADE_DATASET_GET } from '../../modules/local/nextclade_dataset_get'
include { NEXTCLADE_RUN         } from '../../modules/local/nextclade_run'
include { PARSE_NEXTCLADE       } from '../../modules/local/parse_nextclade'

workflow GENOTYPING {
    take:
    consensus       // [meta, fasta]
    serotype_json   // [meta, serotype.json]
    assignment      // [meta, serotype, ref_file, dataset_name]
    coverage_json   // [meta, coverage.json]
    no_file

    main:
    ch_versions = channel.empty()

    ch_assigned_cons = consensus
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(assignment.map { meta, sero, ref_file, ds -> [meta.id, sero, ref_file, ds] },
              failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, fasta, sero, ref_file, ds -> [meta, fasta, sero, ref_file, ds] }

    ch_runnable = ch_assigned_cons
        .filter { _meta, _fasta, _sero, ref_file, ds -> ref_file != 'UNKNOWN' && ds != 'unknown' }
        .map { meta, fasta, _sero, _ref_file, ds -> [ds, meta, fasta] }

    ch_needed_datasets = ch_runnable
        .map { ds, _meta, _fasta -> ds }
        .unique()

    if (params.nextclade_dir) {
        ch_datasets = ch_needed_datasets.map { ds ->
            def short_name = ds.tokenize('/').last()
            def dir = file("${params.nextclade_dir}/${short_name}", checkIfExists: true)
            [ds, dir]
        }
    } else {
        NEXTCLADE_DATASET_GET(ch_needed_datasets)
        ch_datasets = NEXTCLADE_DATASET_GET.out.dataset
        ch_versions = ch_versions.mix(NEXTCLADE_DATASET_GET.out.versions.first())
    }

    ch_nextclade_in = ch_runnable
        .combine(ch_datasets, by: 0)
        .map { _ds, meta, fasta, ds_dir -> [meta, fasta, ds_dir] }

    NEXTCLADE_RUN(ch_nextclade_in)
    ch_versions = ch_versions.mix(NEXTCLADE_RUN.out.versions.first())

    ch_parse_in = consensus
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(serotype_json.map { meta, j -> [meta.id, j] }, failOnMismatch: true, failOnDuplicate: true)
        .join(coverage_json.map { meta, j -> [meta.id, j] }, failOnMismatch: true, failOnDuplicate: true)
        .join(NEXTCLADE_RUN.out.csv.map { meta, csv -> [meta.id, csv] }, remainder: true)
        .map { _id, meta, fasta, sero_j, cov_j, csv -> [meta, fasta, sero_j, cov_j, csv ?: no_file] }

    PARSE_NEXTCLADE(ch_parse_in)
    ch_versions = ch_versions.mix(PARSE_NEXTCLADE.out.versions.first())

    emit:
    json     = PARSE_NEXTCLADE.out.json
    csv      = NEXTCLADE_RUN.out.csv
    versions = ch_versions
}
