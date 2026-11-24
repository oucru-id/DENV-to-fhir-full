process CALL_VARIANTS_ILLUMINA {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::bcftools=1.19"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

    input:
    tuple val(meta), path(bam), path(bai), path(reference)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), path(reference), emit: vcf
    path "versions.yml",                                                                                          emit: versions

    script:
    """
    bcftools mpileup \\
        -a FORMAT/AD,FORMAT/DP \\
        -d 1000000 \\
        -Q 20 \\
        -q 20 \\
        -f ${reference} \\
        ${bam} \\
    | bcftools call --ploidy 1 -mv -Ou \\
    | bcftools norm -f ${reference} -m -any -Ou \\
    | bcftools filter \\
        -i "FORMAT/DP>=${params.consensus_min_depth} && QUAL>=${params.variant_min_quality} && (FORMAT/AD[0:1]/FORMAT/DP)>=${params.consensus_min_af}" \\
        -Oz -o ${meta.id}.filtered.vcf.gz

    bcftools index -t ${meta.id}.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' | gzip > ${meta.id}.filtered.vcf.gz
    touch ${meta.id}.filtered.vcf.gz.tbi
    printf '"%s":\\n    bcftools: stub\\n' "${task.process}" > versions.yml
    """
}
