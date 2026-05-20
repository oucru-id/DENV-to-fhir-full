#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_FHIR {
    publishDir "${params.results_dir}/fhir", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus), path(coverage)
    path(lineage_files)
    each path(org_metadata)

    output:
    path "*.fhir.json", emit: fhir_output
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p lineage_data
    
    for file in ${lineage_files}; do
        if [ -f "\$file" ]; then
            cp "\$file" lineage_data/
        fi
    done

    python3 $baseDir/scripts/annotated_to_fhir.py \\
        --input ${consensus} \\
        --output ${sample_id}.fhir.json \\
        --lineage_dir lineage_data/ \\
        --coverage_file ${coverage} \\
        --organization_metadata ${org_metadata}

    cat <<-END_VERSIONS > versions.yml
    "fhir_converter":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.fhir.json
    touch versions.yml
    """
}

workflow FHIR {
    take:
    fhir_input_ch
    lineage_ch
    org_ch

    main:
    CREATE_FHIR(fhir_input_ch, lineage_ch, org_ch)

    emit:
    fhir_output = CREATE_FHIR.out.fhir_output
    versions    = CREATE_FHIR.out.versions
}
