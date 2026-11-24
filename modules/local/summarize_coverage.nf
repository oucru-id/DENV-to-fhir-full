process SUMMARIZE_COVERAGE {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(per_base_bed), path(summary_txt), path(reference)

    output:
    tuple val(meta), path("${meta.id}.coverage.json"), emit: json
    tuple val(meta), path("${meta.id}.mask.bed"),      emit: mask
    path "versions.yml",                               emit: versions

    script:
    """
    summarize_coverage.py \\
        --mode mosdepth \\
        --sample_id '${meta.id}' \\
        --per_base_bed ${per_base_bed} \\
        --reference ${reference} \\
        --min_depth ${params.consensus_min_depth} \\
        --output_json ${meta.id}.coverage.json \\
        --output_mask_bed ${meta.id}.mask.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id":"${meta.id}","assessed":true,"reference":"NC_001477.1","reference_length":10735,"mean_depth":50.0,"median_depth":50.0,"breadth_1x":1.0,"breadth_min_depth":1.0,"breadth_30x":1.0,"min_depth_threshold":10,"masked_bases":0}' > ${meta.id}.coverage.json
    : > ${meta.id}.mask.bed
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
