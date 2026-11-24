include { BUILD_SEROTYPE_PANEL } from '../../modules/local/build_serotype_panel'

workflow PREPARE_REFERENCE {
    take:
    reference_dir
    reference_panel

    main:
    BUILD_SEROTYPE_PANEL(reference_dir, reference_panel)

    emit:
    panel_index  = BUILD_SEROTYPE_PANEL.out.index
    tagged_fasta = BUILD_SEROTYPE_PANEL.out.tagged_fasta
    versions     = BUILD_SEROTYPE_PANEL.out.versions
}
