include { FASTQC        } from '../../modules/local/fastqc'
include { TRIMMOMATIC   } from '../../modules/local/trimmomatic'
include { NANOPLOT      } from '../../modules/local/nanoplot'
include { CHOPPER       } from '../../modules/local/chopper'
include { HOSTILE_CLEAN } from '../../modules/local/hostile_clean'
include { asFlag        } from './utils_params'

workflow READ_QC {
    take:
    illumina_reads
    ont_reads
    hostile_index

    main:
    ch_versions = channel.empty()

    FASTQC(illumina_reads)
    TRIMMOMATIC(illumina_reads)
    NANOPLOT(ont_reads)
    CHOPPER(ont_reads)

    ch_versions = ch_versions.mix(
        FASTQC.out.versions,
        TRIMMOMATIC.out.versions,
        NANOPLOT.out.versions,
        CHOPPER.out.versions
    )

    if (!asFlag(params.skip_host_removal)) {
        HOSTILE_CLEAN(TRIMMOMATIC.out.reads.mix(CHOPPER.out.reads), hostile_index)
        ch_ill_clean = HOSTILE_CLEAN.out.reads.filter { meta, _r -> meta.platform == 'illumina' }
        ch_ont_clean = HOSTILE_CLEAN.out.reads.filter { meta, _r -> meta.platform == 'ont' }
        ch_versions  = ch_versions.mix(HOSTILE_CLEAN.out.versions)
    } else {
        ch_ill_clean = TRIMMOMATIC.out.reads
        ch_ont_clean = CHOPPER.out.reads
    }

    ch_qc = FASTQC.out.zip.mix(NANOPLOT.out.stats)

    emit:
    illumina_reads = ch_ill_clean
    ont_reads      = ch_ont_clean
    qc             = ch_qc
    versions       = ch_versions
}
