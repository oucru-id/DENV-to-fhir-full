#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TRIMMOMATIC {
    publishDir "${params.results_dir}/qc/trimming", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_{1,2}.fastq.gz"), emit: trimmed_reads
    path "*.log", emit: log

    script:
    """
    java -jar /usr/share/java/trimmomatic.jar PE -threads ${task.cpus} \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_trimmed_1.fastq.gz ${sample_id}_unpaired_1.fastq.gz \\
        ${sample_id}_trimmed_2.fastq.gz ${sample_id}_unpaired_2.fastq.gz \\
        SLIDINGWINDOW:4:20 MINLEN:40 \\
        2> ${sample_id}.trimmomatic.log
    """
}

process CHOPPER {
    publishDir "${params.results_dir}/qc/trimming", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads
    path "*.log", emit: log

    script:
    """
    echo "Running Chopper on ${sample_id}" > ${sample_id}.chopper.log
    
    zcat ${reads} | chopper -q 7 -l 50 | gzip > ${sample_id}_trimmed.fastq.gz 2>> ${sample_id}.chopper.log
    """
}
