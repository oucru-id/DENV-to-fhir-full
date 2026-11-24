include { samplesheetToList  } from 'plugin/nf-schema'
include { PREPARE_REFERENCE  } from '../subworkflows/local/prepare_reference'
include { READ_QC            } from '../subworkflows/local/read_qc'
include { SEROTYPING         } from '../subworkflows/local/serotyping'
include { CONSENSUS_ILLUMINA } from '../subworkflows/local/consensus_illumina'
include { CONSENSUS_NANOPORE } from '../subworkflows/local/consensus_nanopore'
include { GENOTYPING         } from '../subworkflows/local/genotyping'
include { FHIR               } from '../subworkflows/local/fhir'
include { PREPARE_CONSENSUS  } from '../modules/local/prepare_consensus'
include { CONSENSUS_UNMAPPED } from '../modules/local/consensus_unmapped'
include { MULTIQC            } from '../modules/local/multiqc'
include { SAMPLE_REPORT      } from '../modules/local/sample_report'
include { asFlag             } from '../subworkflows/local/utils_params'

def fhirId(String id) {
    return id.replaceAll(/[^A-Za-z0-9\-]/, '-').replaceAll(/-+/, '-').replaceAll(/^-+|-+$/, '')
}

def present(value) {
    return value != null && value != [] && value.toString() != ''
}

def samplesFromSheet(String sheet) {
    def samples = []
    def rows = samplesheetToList(sheet, "${projectDir}/assets/schema_input.json")
    rows.eachWithIndex { row, index ->
        def meta = row[0]
        def fastq_1 = row[1]
        def fastq_2 = row[2]
        def fasta = row[3]
        def line = index + 2
        if (meta.platform == 'illumina') {
            if (!present(fastq_1) || !present(fastq_2)) {
                error "samplesheet line ${line} (${meta.id}): illumina needs fastq_1 and fastq_2"
            }
            samples << [meta: [id: meta.id, platform: 'illumina'], files: [fastq_1, fastq_2]]
        } else if (meta.platform == 'ont') {
            if (!present(fastq_1)) {
                error "samplesheet line ${line} (${meta.id}): ont needs fastq_1"
            }
            samples << [meta: [id: meta.id, platform: 'ont'], files: fastq_1]
        } else {
            if (!present(fasta)) {
                error "samplesheet line ${line} (${meta.id}): fasta needs a fasta column"
            }
            samples << [meta: [id: meta.id, platform: 'fasta'], files: fasta]
        }
    }
    return samples
}

def samplesFromDirectories() {
    def samples = []
    if (!params.reads_dir) {
        return samples
    }
    def pairs = [:]
    files("${params.reads_dir}/*_{1,2}_illumina.fastq.gz").each { f ->
        def id = f.name.replaceFirst(/_[12]_illumina\.fastq\.gz$/, '')
        if (!pairs.containsKey(id)) {
            pairs[id] = []
        }
        pairs[id] << f
    }
    pairs.keySet().sort().each { id ->
        def reads = pairs[id].sort { f -> f.name }
        if (reads.size() != 2) {
            error "Illumina sample ${id}: expected _1_illumina and _2_illumina files, found ${reads*.name}"
        }
        samples << [meta: [id: id, platform: 'illumina'], files: reads]
    }
    files("${params.reads_dir}/*_ont.fastq.gz").sort { f -> f.name }.each { f ->
        samples << [meta: [id: f.name.replaceFirst(/_ont\.fastq\.gz$/, ''), platform: 'ont'], files: f]
    }
    def fastas = files("${params.reads_dir}/*.fasta") + files("${params.reads_dir}/*.fa") + files("${params.reads_dir}/*.fna")
    fastas.sort { f -> f.name }.each { f ->
        samples << [meta: [id: f.name.replaceFirst(/\.(fasta|fa|fna)$/, ''), platform: 'fasta'], files: f]
    }
    return samples
}

def validateSamples(samples) {
    if (!samples) {
        error "No samples found. Give --input samplesheet.csv, or --reads_dir " +
              "(files named <id>_{1,2}_illumina.fastq.gz, <id>_ont.fastq.gz or <id>.fasta)."
    }
    def seen = [:]
    def seenFhir = [:]
    samples.each { s ->
        def id = s.meta.id
        if (!(id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]{0,39}$/)) {
            error "Sample id '${id}': use 1-40 letters, digits, '.', '_' or '-', starting with a letter or digit"
        }
        if (seen.containsKey(id)) {
            error "Sample id '${id}' appears more than once"
        }
        def fid = fhirId(id)
        if (seenFhir.containsKey(fid)) {
            error "Sample ids '${seenFhir[fid]}' and '${id}' become the same FHIR id '${fid}'; rename one"
        }
        if (id ==~ /^\d{16}$/) {
            log.warn "Sample id '${id}' looks like a 16-digit national ID (NIK). Sample ids appear in " +
                     "the PHI-free genomic output; use a laboratory accession instead."
        }
        seen[id] = true
        seenFhir[fid] = id
    }
    if (samples.any { s -> s.meta.platform == 'ont' } && !params.medaka_model) {
        error "ONT samples need --medaka_model (e.g. r1041_e82_400bps_sup_variant_v5.0.0); " +
              "list the available models with `medaka tools list_models`."
    }
    if (params.fhir_server_url && !params.patient_metadata) {
        error "--fhir_server_url uploads linked bundles, which need --patient_metadata"
    }
}

workflow DENVFHIR {
    main:
    def samples = params.input ? samplesFromSheet(params.input) : samplesFromDirectories()
    validateSamples(samples)
    def counts = samples.countBy { s -> s.meta.platform }
    log.info "Samples: ${samples.size()} (illumina: ${counts.illumina ?: 0}, ont: ${counts.ont ?: 0}, fasta: ${counts.fasta ?: 0})"

    def no_file = file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    def reference_dir = file(params.reference_dir, checkIfExists: true)
    def reference_panel = file(params.reference_panel, checkIfExists: true)
    def hostile_index = params.hostile_index ? file(params.hostile_index, checkIfExists: true) : no_file
    def primer_bed = params.primer_bed ? file(params.primer_bed, checkIfExists: true) : no_file

    ch_samples  = channel.fromList(samples).map { s -> [s.meta, s.files] }
    ch_illumina = ch_samples.filter { meta, _f -> meta.platform == 'illumina' }
    ch_ont      = ch_samples.filter { meta, _f -> meta.platform == 'ont' }
    ch_fasta    = ch_samples.filter { meta, _f -> meta.platform == 'fasta' }
    ch_versions = channel.empty()

    PREPARE_REFERENCE(reference_dir, reference_panel)
    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    READ_QC(ch_illumina, ch_ont, hostile_index)
    PREPARE_CONSENSUS(ch_fasta)
    ch_versions = ch_versions.mix(READ_QC.out.versions, PREPARE_CONSENSUS.out.versions)

    ch_for_serotype = READ_QC.out.illumina_reads
        .mix(READ_QC.out.ont_reads, PREPARE_CONSENSUS.out.consensus)
    SEROTYPING(ch_for_serotype, PREPARE_REFERENCE.out.panel_index, reference_panel)
    ch_versions = ch_versions.mix(SEROTYPING.out.versions)

    ch_assign_ill = SEROTYPING.out.assignment
        .filter { meta, _s, _r, _d -> meta.platform == 'illumina' }
        .map { meta, serotype, ref_file, dataset -> [meta.id, serotype, ref_file, dataset] }

    ch_assign_ont = SEROTYPING.out.assignment
        .filter { meta, _s, _r, _d -> meta.platform == 'ont' }
        .map { meta, serotype, ref_file, dataset -> [meta.id, serotype, ref_file, dataset] }

    ch_ill_joined = READ_QC.out.illumina_reads
        .map { meta, reads -> [meta.id, meta, reads] }
        .join(ch_assign_ill, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, reads, serotype, ref_file, _ds -> [meta, reads, serotype, ref_file] }

    ch_ont_joined = READ_QC.out.ont_reads
        .map { meta, reads -> [meta.id, meta, reads] }
        .join(ch_assign_ont, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, reads, serotype, ref_file, _ds -> [meta, reads, serotype, ref_file] }

    ch_ill_mapped = ch_ill_joined
        .filter { _meta, _reads, _sero, ref_file -> ref_file != 'UNKNOWN' }
        .map { meta, reads, serotype, ref_file ->
            [meta, reads, file("${params.reference_dir}/${ref_file}", checkIfExists: true), serotype]
        }

    ch_ont_mapped = ch_ont_joined
        .filter { _meta, _reads, _sero, ref_file -> ref_file != 'UNKNOWN' }
        .map { meta, reads, serotype, ref_file ->
            [meta, reads, file("${params.reference_dir}/${ref_file}", checkIfExists: true), serotype]
        }

    ch_unmapped_meta = ch_ill_joined.mix(ch_ont_joined)
        .filter { _meta, _reads, _sero, ref_file -> ref_file == 'UNKNOWN' }
        .map { meta, _reads, _sero, _ref -> meta }

    CONSENSUS_ILLUMINA(ch_ill_mapped, primer_bed)
    CONSENSUS_NANOPORE(ch_ont_mapped, primer_bed)
    CONSENSUS_UNMAPPED(ch_unmapped_meta)
    ch_versions = ch_versions.mix(
        CONSENSUS_ILLUMINA.out.versions,
        CONSENSUS_NANOPORE.out.versions,
        CONSENSUS_UNMAPPED.out.versions
    )

    ch_consensus = CONSENSUS_ILLUMINA.out.consensus
        .mix(CONSENSUS_NANOPORE.out.consensus, PREPARE_CONSENSUS.out.consensus, CONSENSUS_UNMAPPED.out.consensus)
    ch_coverage = CONSENSUS_ILLUMINA.out.coverage
        .mix(CONSENSUS_NANOPORE.out.coverage, PREPARE_CONSENSUS.out.coverage, CONSENSUS_UNMAPPED.out.coverage)

    GENOTYPING(ch_consensus, SEROTYPING.out.json, SEROTYPING.out.assignment, ch_coverage, no_file)
    ch_versions = ch_versions.mix(GENOTYPING.out.versions)

    ch_qc = READ_QC.out.qc
        .mix(CONSENSUS_ILLUMINA.out.mqc, CONSENSUS_NANOPORE.out.mqc)
        .flatMap { _meta, qc_files -> [qc_files].flatten() }
    MULTIQC(ch_qc.collect().filter { list -> !list.isEmpty() })
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    ch_report_in = ch_consensus
        .map { meta, fasta -> [meta.id, meta, fasta] }
        .join(SEROTYPING.out.json.map { meta, j -> [meta.id, j] }, failOnMismatch: true, failOnDuplicate: true)
        .join(GENOTYPING.out.json.map { meta, j -> [meta.id, j] }, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_coverage.map { meta, j -> [meta.id, j] }, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, fasta, sero_j, geno_j, cov_j -> [meta, fasta, sero_j, geno_j, cov_j] }

    SAMPLE_REPORT(ch_report_in)
    ch_versions = ch_versions.mix(SAMPLE_REPORT.out.versions)

    FHIR(ch_consensus, SEROTYPING.out.json, GENOTYPING.out.json, ch_coverage, no_file)
    ch_versions = ch_versions.mix(FHIR.out.versions)

    def versions_header = ("\"Workflow\":\n    DENVtoFHIR: ${workflow.manifest.version}\n" +
                           "    Nextflow: ${nextflow.version}\n" +
                           "    fhir_profile: ${params.fhir_profile}\n" +
                           "    primer_trimmed: ${params.primer_bed ? 'true' : 'false'}\n" +
                           "    skip_host_removal: ${asFlag(params.skip_host_removal)}\n" +
                           "    reference_panel: ${reference_panel.name}\n").toString()
    ch_versions
        .unique()
        .collectFile(name: 'versions_collected.yml', sort: true, newLine: true)
        .map { collected -> versions_header + collected.text }
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")

    emit:
    consensus = ch_consensus
    genomic   = FHIR.out.genomic
    bundles   = FHIR.out.bundles
}
