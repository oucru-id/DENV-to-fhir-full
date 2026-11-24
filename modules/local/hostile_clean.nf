process HOSTILE_CLEAN {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::hostile=2.0.2"
    container "quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0"

    input:
    tuple val(meta), path(reads)
    path hostile_index

    output:
    tuple val(meta), path("${meta.id}*.clean.fastq.gz"), emit: reads
    path "${meta.id}.hostile.log",                       emit: log
    path "versions.yml",                                 emit: versions

    script:
    def cache_export = hostile_index.name != 'NO_FILE' ? "export HOSTILE_CACHE_DIR=\"\$PWD/${hostile_index}\"" : ""
    if (meta.platform == 'illumina') {
        """
        ${cache_export}
        hostile clean \\
            --fastq1 ${reads[0]} \\
            --fastq2 ${reads[1]} \\
            --aligner bowtie2 \\
            --index ${params.hostile_index_name} \\
            --threads ${task.cpus} \\
            --output hostile_out \\
            > ${meta.id}.hostile.log 2>&1

        f1=\$(find hostile_out -maxdepth 1 -name '*clean_1*.fastq.gz' -o -name '*_1.clean*.fastq.gz' | head -n 1)
        f2=\$(find hostile_out -maxdepth 1 -name '*clean_2*.fastq.gz' -o -name '*_2.clean*.fastq.gz' | head -n 1)
        if [ -z "\$f1" ] || [ -z "\$f2" ]; then
            echo "ERROR: Hostile did not produce paired clean FASTQs" >&2
            ls -la hostile_out >&2
            exit 1
        fi
        mv "\$f1" "${meta.id}_1.clean.fastq.gz"
        mv "\$f2" "${meta.id}_2.clean.fastq.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version 2>&1 | tail -n 1)
        END_VERSIONS
        """
    } else {
        """
        ${cache_export}
        hostile clean \\
            --fastq1 ${reads} \\
            --aligner minimap2 \\
            --index ${params.hostile_index_name} \\
            --threads ${task.cpus} \\
            --output hostile_out \\
            > ${meta.id}.hostile.log 2>&1

        f1=\$(find hostile_out -maxdepth 1 -name '*clean*.fastq.gz' | head -n 1)
        if [ -z "\$f1" ]; then
            echo "ERROR: Hostile did not produce clean FASTQ" >&2
            ls -la hostile_out >&2
            exit 1
        fi
        mv "\$f1" "${meta.id}.clean.fastq.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hostile: \$(hostile --version 2>&1 | tail -n 1)
        END_VERSIONS
        """
    }

    stub:
    if (meta.platform == 'illumina') {
        """
        echo | gzip > ${meta.id}_1.clean.fastq.gz
        echo | gzip > ${meta.id}_2.clean.fastq.gz
        touch ${meta.id}.hostile.log
        printf '"%s":\\n    hostile: stub\\n' "${task.process}" > versions.yml
        """
    } else {
        """
        echo | gzip > ${meta.id}.clean.fastq.gz
        touch ${meta.id}.hostile.log
        printf '"%s":\\n    hostile: stub\\n' "${task.process}" > versions.yml
        """
    }
}
