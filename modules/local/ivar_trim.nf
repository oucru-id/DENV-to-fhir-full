include { asFlag } from '../../subworkflows/local/utils_params'

process IVAR_TRIM {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::ivar=1.4.4 bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/ivar:1.4.4--h077b44d_0"

    input:
    tuple val(meta), path(bam), path(bai), path(reference), val(serotype)
    path primer_bed

    output:
    tuple val(meta), path("${meta.id}.trimmed.bam"), path("${meta.id}.trimmed.bam.bai"), path(reference), emit: bam
    path "versions.yml",                                                                                  emit: versions

    script:
    def include_flag = asFlag(params.ivar_include_no_primer) ? '-e' : ''
    def sero_nodash  = serotype.replace('-', '')
    """
    BED=""
    if [ -d "${primer_bed}" ]; then
        for cand in "${primer_bed}/${serotype}.bed" "${primer_bed}/${sero_nodash}.bed" "${primer_bed}/${serotype.toLowerCase()}.bed" "${primer_bed}/${sero_nodash.toLowerCase()}.bed"; do
            if [ -f "\$cand" ]; then
                BED="\$cand"
                break
            fi
        done
        if [ -z "\$BED" ]; then
            echo "ERROR: no primer BED for serotype ${serotype} in ${primer_bed}" >&2
            exit 1
        fi
    else
        BED="${primer_bed}"
    fi

    ivar trim \\
        -b "\$BED" \\
        -p ${meta.id}.ivar \\
        -i ${bam} \\
        -q ${params.ivar_min_quality} \\
        -m ${params.ivar_min_length} \\
        ${include_flag}

    samtools sort -@ ${task.cpus} -o ${meta.id}.trimmed.bam ${meta.id}.ivar.bam
    samtools index ${meta.id}.trimmed.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version | head -n 1 | sed 's/iVar version //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.trimmed.bam ${meta.id}.trimmed.bam.bai
    printf '"%s":\\n    ivar: stub\\n' "${task.process}" > versions.yml
    """
}
