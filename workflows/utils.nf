nextflow.enable.dsl = 2

process VERSIONS {
    publishDir "${params.results_dir}", mode: 'copy'

    output:
    path "software_versions.yml"

    script:
    """
    echo "  name: dengue_pipeline v.${params.version}" >> software_versions.yml
    echo "  nextflow: $nextflow.version" >> software_versions.yml
    
    echo "databases:" >> software_versions.yml
    echo "  reference_dir: ${params.reference_dir}" >> software_versions.yml
    
    echo "processing_settings:" >> software_versions.yml
    echo "  illumina:" >> software_versions.yml
    echo "    trimmomatic: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36'" >> software_versions.yml
    echo "  nanopore:" >> software_versions.yml
    echo "    chopper_min_q: 7" >> software_versions.yml
    echo "    chopper_min_l: 50" >> software_versions.yml
    echo "  genotyping:" >> software_versions.yml
    
    export BASE_DIR="$baseDir"
    python3 $baseDir/scripts/get_versions.py >> software_versions.yml
    """
}