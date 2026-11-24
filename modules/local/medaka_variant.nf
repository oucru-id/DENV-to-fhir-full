process MEDAKA_VARIANT {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::medaka=2.1.0"
    container "quay.io/biocontainers/medaka:2.1.0--py38ha0c3a46_0"

    input:
    tuple val(meta), path(bam), path(bai), path(reference)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), path(reference), emit: vcf
    path "versions.yml",                                                                                          emit: versions

    script:
    """
    cp -L ${reference} reference.fasta
    samtools faidx reference.fasta
    export CUDA_VISIBLE_DEVICES=""

    medaka inference ${bam} consensus_probs.hdf \\
        --model ${params.medaka_model} \\
        --batch_size ${params.medaka_batch_size} \\
        --threads ${task.cpus}
    medaka vcf consensus_probs.hdf reference.fasta medaka.vcf
    bcftools sort medaka.vcf -o medaka.sorted.vcf
    medaka tools annotate medaka.sorted.vcf reference.fasta ${bam} medaka.annotated.vcf

    if [ ! -s medaka.annotated.vcf ]; then
        echo "ERROR: medaka produced no annotated VCF" >&2
        exit 1
    fi

    bcftools norm -f reference.fasta -m -any medaka.annotated.vcf -Ou \\
    | bcftools view -v snps,indels -i 'GT="1" || GT="1/1" || GT="0/1"' -Ou \\
    | bcftools filter \\
        -i "INFO/DP>=${params.consensus_min_depth} && (QUAL>=${params.variant_min_quality} || FORMAT/GQ>=${params.variant_min_quality})" \\
        -Oz -o ${meta.id}.filtered.vcf.gz

    bcftools index -t ${meta.id}.filtered.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version 2>&1 | sed 's/medaka //')
        medaka_model: ${params.medaka_model}
    END_VERSIONS
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' | gzip > ${meta.id}.filtered.vcf.gz
    touch ${meta.id}.filtered.vcf.gz.tbi
    printf '"%s":\\n    medaka: stub\\n' "${task.process}" > versions.yml
    """
}
