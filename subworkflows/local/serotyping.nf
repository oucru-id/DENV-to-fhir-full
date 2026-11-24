include { SEROTYPE_MAP      } from '../../modules/local/serotype_map'
include { SEROTYPE_CLASSIFY } from '../../modules/local/serotype_classify'

workflow SEROTYPING {
    take:
    reads_or_fasta
    panel_index
    reference_panel

    main:
    SEROTYPE_MAP(reads_or_fasta, panel_index)
    SEROTYPE_CLASSIFY(SEROTYPE_MAP.out.paf, reference_panel)

    emit:
    json       = SEROTYPE_CLASSIFY.out.json
    assignment = SEROTYPE_CLASSIFY.out.assignment
    versions   = SEROTYPE_MAP.out.versions.mix(SEROTYPE_CLASSIFY.out.versions)
}
