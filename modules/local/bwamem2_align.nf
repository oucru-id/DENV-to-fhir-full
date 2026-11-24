process BWAMEM2_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0"

    input:
    tuple val(meta), path(reads), path(reference)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), path(reference), emit: bam
    path "versions.yml",                                                                  emit: versions

    script:
    """
    cp -L ${reference} ref.fasta
    bwa-mem2 index ref.fasta

    bwa-mem2 mem -t ${task.cpus} \\
        -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:lib1" \\
        ref.fasta ${reads[0]} ${reads[1]} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null | tail -1)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    printf '"%s":\\n    bwa-mem2: stub\\n' "${task.process}" > versions.yml
    """
}
