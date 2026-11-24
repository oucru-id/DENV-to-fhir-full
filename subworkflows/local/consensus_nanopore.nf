include { MINIMAP2_ALIGN     } from '../../modules/local/minimap2_align'
include { IVAR_TRIM          } from '../../modules/local/ivar_trim'
include { MOSDEPTH           } from '../../modules/local/mosdepth'
include { MEDAKA_VARIANT     } from '../../modules/local/medaka_variant'
include { SUMMARIZE_COVERAGE } from '../../modules/local/summarize_coverage'
include { BUILD_CONSENSUS    } from '../../modules/local/build_consensus'

workflow CONSENSUS_NANOPORE {
    take:
    assigned   
    primer_bed

    main:
    ch_versions = channel.empty()

    MINIMAP2_ALIGN(assigned.map { meta, reads, ref, _sero -> [meta, reads, ref] })
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    if (primer_bed.name != 'NO_FILE') {
        ch_trim_in = MINIMAP2_ALIGN.out.bam
            .map { meta, bam, bai, ref -> [meta.id, meta, bam, bai, ref] }
            .join(assigned.map { meta, _r, _ref, sero -> [meta.id, sero] }, failOnMismatch: true, failOnDuplicate: true)
            .map { _id, meta, bam, bai, ref, sero -> [meta, bam, bai, ref, sero] }
        IVAR_TRIM(ch_trim_in, primer_bed)
        ch_bam = IVAR_TRIM.out.bam
        ch_versions = ch_versions.mix(IVAR_TRIM.out.versions)
    } else {
        ch_bam = MINIMAP2_ALIGN.out.bam
    }

    MOSDEPTH(ch_bam)
    MEDAKA_VARIANT(ch_bam)
    SUMMARIZE_COVERAGE(MOSDEPTH.out.coverage)

    ch_cons_in = MEDAKA_VARIANT.out.vcf
        .map { meta, vcf, tbi, ref -> [meta.id, meta, vcf, tbi, ref] }
        .join(SUMMARIZE_COVERAGE.out.mask.map { meta, mask -> [meta.id, mask] }, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, vcf, tbi, ref, mask -> [meta, vcf, tbi, ref, mask] }

    BUILD_CONSENSUS(ch_cons_in)

    ch_versions = ch_versions.mix(
        MOSDEPTH.out.versions,
        MEDAKA_VARIANT.out.versions,
        SUMMARIZE_COVERAGE.out.versions,
        BUILD_CONSENSUS.out.versions
    )

    emit:
    consensus = BUILD_CONSENSUS.out.consensus
    vcf       = MEDAKA_VARIANT.out.vcf
    bam       = ch_bam
    coverage  = SUMMARIZE_COVERAGE.out.json
    mqc       = MOSDEPTH.out.mqc
    versions  = ch_versions
}
