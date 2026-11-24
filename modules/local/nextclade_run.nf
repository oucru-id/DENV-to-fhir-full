process NEXTCLADE_RUN {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::nextclade=3.29.1"
    container "quay.io/biocontainers/nextclade:3.29.1--h9ee0642_0"

    input:
    tuple val(meta), path(consensus), path(dataset_dir)

    output:
    tuple val(meta), path("${meta.id}.nextclade.csv"),  emit: csv
    tuple val(meta), path("${meta.id}.nextclade.json"), emit: json
    path "versions.yml",                                emit: versions

    script:
    """
    nextclade run \\
        --jobs ${task.cpus} \\
        --input-dataset ${dataset_dir} \\
        --output-csv ${meta.id}.nextclade.csv \\
        --output-json ${meta.id}.nextclade.json \\
        ${consensus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(nextclade --version | sed 's/nextclade //')
    END_VERSIONS
    """

    stub:
    """
    printf 'index;seqName;clade;qc.overallScore;qc.overallStatus;coverage;aaSubstitutions;aaDeletions;aaInsertions;substitutions;deletions;insertions;frameShifts;totalMissing;totalNonACGTNs\\n0;${meta.id};1V_D.1.1;0;good;0.99;E:D37N;;;A100G;;;;0;0\\n' > ${meta.id}.nextclade.csv
    echo '{}' > ${meta.id}.nextclade.json
    printf '"%s":\\n    nextclade: stub\\n' "${task.process}" > versions.yml
    """
}
