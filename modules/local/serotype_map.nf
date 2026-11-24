process SEROTYPE_MAP {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::minimap2=2.26"
    container "quay.io/biocontainers/minimap2:2.26--he4a0461_2"

    input:
    tuple val(meta), path(reads)
    path panel_index

    output:
    tuple val(meta), path("${meta.id}.minimap2.paf"), env('TOTAL_INPUT_READS'), emit: paf
    path "versions.yml",                                                       emit: versions

    script:
    def is_illumina = (meta.platform == 'illumina')
    def is_fasta    = (meta.platform == 'fasta')
    def read_file   = is_illumina ? reads[0] : reads
    def mm2_preset  = is_fasta ? 'asm20' : (is_illumina ? 'sr' : 'map-ont')
    if (is_fasta) {
        """
        TOTAL_INPUT_READS=\$(grep -c '^>' ${read_file} || echo 0)
        if [ "\$TOTAL_INPUT_READS" -gt 0 ]; then
            minimap2 -x ${mm2_preset} --secondary=no -t ${task.cpus} \\
                ${panel_index} ${read_file} > ${meta.id}.minimap2.paf || touch ${meta.id}.minimap2.paf
        else
            : > ${meta.id}.minimap2.paf
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version)
        END_VERSIONS
        """
    } else {
        """
        MAX_LINES=\$((180000 * 4))
        if [[ "${read_file}" == *.gz ]]; then
            (zcat "${read_file}" 2>/dev/null || true) | head -n \$MAX_LINES > subset.fastq
        else
            head -n \$MAX_LINES "${read_file}" > subset.fastq
        fi

        TOTAL_INPUT_READS=\$(( \$(wc -l < subset.fastq) / 4 ))
        if [ "\$TOTAL_INPUT_READS" -gt 0 ]; then
            minimap2 -x ${mm2_preset} --secondary=no -t ${task.cpus} \\
                ${panel_index} subset.fastq > ${meta.id}.minimap2.paf || touch ${meta.id}.minimap2.paf
        else
            : > ${meta.id}.minimap2.paf
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version)
        END_VERSIONS
        """
    }

    stub:
    """
    touch ${meta.id}.minimap2.paf
    TOTAL_INPUT_READS=1000
    printf '"%s":\\n    minimap2: stub\\n' "${task.process}" > versions.yml
    """
}
