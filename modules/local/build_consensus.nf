process BUILD_CONSENSUS {
    tag "${meta.id}"
    label 'process_low'
    conda "bioconda::bcftools=1.19"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

    input:
    tuple val(meta), path(vcf), path(tbi), path(reference), path(mask_bed)

    output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: consensus
    path "versions.yml",                                 emit: versions

    script:
    """
    MASK_ARG=""
    if [ -s "${mask_bed}" ]; then
        MASK_ARG="-m ${mask_bed}"
    fi

    bcftools consensus \\
        -f ${reference} \\
        \$MASK_ARG \\
        ${vcf} \\
    | sed "1s/^>.*/>${meta.id}/" \\
    > ${meta.id}.consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    printf '>${meta.id}\\nACGTACGTACGT\\n' > ${meta.id}.consensus.fasta
    printf '"%s":\\n    bcftools: stub\\n' "${task.process}" > versions.yml
    """
}
