#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MULTIQC {
    publishDir "${params.results_dir}/qc", mode: 'copy'

    input:
    path(qc_files)

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p fastqc_results
    cp -L ${qc_files} fastqc_results/

    multiqc fastqc_results \
        --force \
        --verbose \
        --outdir . \
        --no-ansi \
        --interactive \
        --zip-data-dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}

process CREATE_DENGUE_REPORTS {
    publishDir "${params.results_dir}/reports", mode: 'copy'
    
    input:
    tuple val(sample_id), path(serotype_json), path(genotype_json), path(consensus)
    
    output:
    path "${sample_id}_dengue_report.txt", emit: dengue_report
    path "versions.yml", emit: versions
    
    script:
    """
    if [ ! -f "${genotype_json}" ] || [ "${genotype_json}" = "empty.json" ]; then
        echo '{}' > empty_genotype.json
        GENOTYPE_ARG=""
    else
        GENOTYPE_ARG="--genotype_json ${genotype_json}"
    fi
    
    python3 $baseDir/scripts/generate_dengue_report.py \
        --sample_id ${sample_id} \
        --serotype_json ${serotype_json} \
        \$GENOTYPE_ARG \
        --consensus ${consensus} \
        --output ${sample_id}_dengue_report.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}

workflow GENERATE_REPORT {
    take:
    qc_files

    main:
    MULTIQC(qc_files)

    emit:
    report = MULTIQC.out.report
    versions = MULTIQC.out.versions
}

workflow GENERATE_DENGUE_REPORTS {
    take:
    combined_data

    main:
    CREATE_DENGUE_REPORTS(combined_data)

    emit:
    reports = CREATE_DENGUE_REPORTS.out.dengue_report
    versions = CREATE_DENGUE_REPORTS.out.versions
}
