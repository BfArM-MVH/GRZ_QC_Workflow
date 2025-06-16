/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQC                         } from '../modules/nf-core/fastqc'
include { FASTP                          } from '../modules/nf-core/fastp'
include { MULTIQC                        } from '../modules/nf-core/multiqc'
include { CONVERT_BED_CHROM              } from '../modules/local/convert_bed_chrom'
include { COMPARE_THRESHOLD              } from '../modules/local/compare_threshold'
include { MERGE_REPORTS                  } from '../modules/local/merge_reports'
include { BWAMEM2_INDEX                  } from '../modules/nf-core/bwamem2/index'
include { MOSDEPTH                       } from '../modules/nf-core/mosdepth'
include { SAMTOOLS_FAIDX                 } from '../modules/nf-core/samtools/faidx'
include { SAVE_REFERENCE                 } from '../modules/local/save_reference'
include { FASTQ_ALIGN_BWA_MARKDUPLICATES } from '../subworkflows/local/fastq_align_bwa_markduplicates'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GRZQC {
    take:
    ch_samplesheet // channel: samplesheet created by parsing metadata.json file
    ch_genome

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // conflicting parameters
    //
    // Error if BWA or FAI is given without FASTA and reference_path is not present
    if ((params.bwa || params.fai) && !params.fasta && !params.reference_path) {
        println("\033[1;31mERROR:\033[0m 'bwa' or 'fai' is provided, but 'fasta' is not present. Please provide a valid FASTA file.")
        System.exit(1)
    }
    // // Warning if reference_path is provided together with fasta, fai, or bwa
    // if ((params.reference_path && params.fasta) || (params.reference_path && params.fai) || (params.reference_path && params.bwa)) {
    //     println("\033[1;33mWARNING:\033[0m 'reference_path' is provided together with 'fasta', 'fai', or 'bwa'. Only 'reference_path' will be considered.")
    // }
    // // Specific warning if reference_path is given with bwa or fai, but fasta is missing
    // if (params.reference_path && (params.bwa || params.fai) && !params.fasta) {
    //     println("\033[1;33mWARNING:\033[0m 'reference_path' is provided together with 'fasta', 'fai', or 'bwa'. Only 'reference_path' will be considered.")
    // }

    //
    // set up channels
    //
    // match fa and fasta extensions
    fastaExts = ['.fa', '.fasta', '.fa.gz', '.fasta.gz']
    faiExts = ['.fa.fai', '.fasta.fai', '.fa.gz.fai', '.fasta.gz.fai']

    // create reference channels
    if (params.reference_path) {

        fasta = ch_genome
            .map { genome ->
                // 1) Look for any FASTA under reference_path/genome
                candidates = fastaExts.collect { ext -> file("${params.reference_path}/${genome}/*${ext}") }.flatten()
                f = candidates.find { it.exists() }
                if (f) {
                    tuple([id: f.baseName], f)
                }
                else {
                    // 2) If user passed --fasta, use that (error if missing)
                    if( params.fasta ) {
                        def f2 = file(params.fasta)
                        if( ! f2.exists() ) { error "User-specified FASTA not found: ${params.fasta}" }
                        return tuple([ id: f2.baseName ], f2)
                    }
                    else {
                        // 3) Otherwise download the default for this genome name
                        defaultFasta = genome == 'GRCh38'
                            ? "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz"
                            : "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
                        def f3 = file(defaultFasta)
                        return tuple([ id: f3.baseName ], f3)
                    }
                }
            }
            .collect()

    }

    //
    // build BWA index and FAI under reference_path
    //
    genomeStr = ch_genome.toString()

    // BWA index
    refBwa = params.reference_path ? file("${params.reference_path}/${genomeStr}/bwamem2") : null

    // DEBUG prints
    println ">>> DEBUG: params.reference_path = ${params.reference_path}"
    println ">>> DEBUG: genomeStr         = ${genomeStr}"
    println ">>> DEBUG: checking refBwa at ${refBwa}"
    // turn it into a File object (or null)
    println ">>> DEBUG: refBwa is $refBwa, exists? ${refBwa?.exists()}"

    if( refBwa && refBwa.exists() ) {
        // 1) existing bwamem2 under reference_path
        println ">>> DEBUG: using existing bwa index"
        bwa = Channel
            .fromPath(refBwa.toString(), checkIfExists: true)
            .map { f -> tuple([id: f.baseName], f) }
    }
    else if( params.bwa ) {
        // 2) user-provided bwa binary
        bwa = Channel
            .fromPath(params.bwa, checkIfExists: true)
            .map { f -> tuple([id: 'bwa'], f) }
    }
    else {
        // 3) fallback: build index via your process
        println ">>> DEBUG: fallback: build index via your process"
        BWAMEM2_INDEX(
            fasta
        )
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
        bwa = BWAMEM2_INDEX.out.index
    }

    // fai
    reffai = params.reference_path ? file("${params.reference_path}/${genomeStr}/*.fai")[0] : null
    if( reffai && reffai.exists() ) {
        // 1) existing fai under reference_path
        fai = Channel
            .fromPath(reffai.toString(), checkIfExists: true)
            .map { f -> tuple([id: f.baseName], f) }
    }
    else if( params.fai ) {
        // 2) user-provided fai
        fai = Channel
            .fromPath(params.fai, checkIfExists: true)
            .map { f -> tuple([id: 'fai'], f) }
    }
    else {
        // 3) fallback: build index via your process
        SAMTOOLS_FAIDX(
            fasta,
            [[], []],
            false,
        )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        fai = SAMTOOLS_FAIDX.out.fai
    }

    // if (!params.reference_path) {

    //     // if (!params.bwa) {
    //     //     // create bwa index if not provided
    //     //     BWAMEM2_INDEX(
    //     //         fasta
    //     //     )

    //     //     ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    //     //     bwa = BWAMEM2_INDEX.out.index
    //     // }
    //     // if (!params.fai) {
    //     //     // create fai index if not provided
    //     //     SAMTOOLS_FAIDX(
    //     //         fasta,
    //     //         [[], []],
    //     //         false,
    //     //     )

    //     //     fai = SAMTOOLS_FAIDX.out.fai
    //     //     ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    //     // }

    //     if (params.save_reference) {
    //         // save reference for the first run
    //         SAVE_REFERENCE(
    //             fasta,
    //             fai,
    //             bwa,
    //             ch_genome,
    //         )
    //     }
    // }

    // TARGET BED channel
    if (params.target) {
        ch_target = Channel.fromPath(params.target, checkIfExists: true).collect()
    }
    else {
        ch_target = ch_genome
            .flatMap { genome ->
                def defaultTargetCh = genome == 'GRCh38'
                    ? "${projectDir}/assets/default_files/hg38_440_omim_genes.bed"
                    : "${projectDir}/assets/default_files/hg19_439_omim_genes.bed"
                def f = file(defaultTargetCh)
                if (!f.exists()) {
                    error("Default target BED missing: ${f}")
                }
                return f
            }
            .collect()
    }

    // CHROMOSOME‐NAME MAPPING channel
    if (params.mapping_chrom) {
        mapping_chrom = Channel.fromPath(params.mapping_chrom, checkIfExists: true).collect()
    }
    else {
        mapping_chrom = ch_genome
            .flatMap { genome ->
                def defaultMappingCh = genome == 'GRCh38'
                    ? "${projectDir}/assets/default_files/hg38_NCBI2UCSC.txt"
                    : "${projectDir}/assets/default_files/hg19_NCBI2UCSC.txt"
                def f = file(defaultMappingCh)
                if (!f.exists()) {
                    error("Default mapping-chrom file missing: ${f}")
                }
                return f
            }
            .collect()
    }

    ch_samplesheet
        .branch { meta, reads, alignment ->
            reads: reads.size() > 0
            return [meta, reads]
            alignments: alignment.size() > 0
            return [meta, alignment]
        }
        .set { samplesheet_ch }

    // Run FASTQC on FASTQ files - per lane
    FASTQC(
        samplesheet_ch.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // Run FASP on FASTQ files - per lane
    save_trimmed_fail = false
    save_merged = false
    FASTP(
        samplesheet_ch.reads,
        [],
        false,
        save_trimmed_fail,
        save_merged,
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { _meta, json -> json })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect { _meta, html -> html })
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // align FASTQs per lane, merge, and sort
    FASTQ_ALIGN_BWA_MARKDUPLICATES(
        samplesheet_ch.reads,
        samplesheet_ch.alignments,
        bwa,
        true,
        fasta,
        fai,
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.stat.collect { _meta, file -> file })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.flagstat.collect { _meta, file -> file })

    ch_bams = FASTQ_ALIGN_BWA_MARKDUPLICATES.out.bam.join(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.bai, by: 0)

    // prepare mosdepth inputs
    // for WGS: defined bed files of ~400 genes
    // for WES and panel: extract bed from samplesheet and
    //      do the conversion with the correct mapping file (if it is NCBI format, covert it to UCSC).
    // for both cases different file required for hg38 and hg19

    ch_bams
        .branch { meta, bam, bai ->
            targeted: meta.libraryType in ["wes", "panel", "wes_lr", "panel_lr"]
            return [meta, bam, bai, meta.bed_file]
            wgs: meta.libraryType in ["wgs", "wgs_lr"]
            return [meta, bam, bai]
        }
        .set { ch_bams_bed }

    // convert given bed files to UCSC style names
    // for WES and panel, run the conversion process: if the bed file has NCBI-style names, they will be converted.
    CONVERT_BED_CHROM(
        ch_bams_bed.targeted.map { meta, _bam, _bai, bed_file -> [meta, bed_file] },
        mapping_chrom,
    )
    ch_converted_bed = CONVERT_BED_CHROM.out.converted_bed
    ch_versions = ch_versions.mix(CONVERT_BED_CHROM.out.versions)

    ch_bams_bed.targeted
        .join(ch_converted_bed)
        .map { meta, bam, bai, _old_bed, converted_bed ->
            def newMeta = meta.clone()
            newMeta.remove('bed_file')
            [newMeta, bam, bai, converted_bed]
        }
        .set { ch_bams_bed_targeted }

    ch_bams_bed.wgs
        .combine(ch_target)
        .map { meta, bam, bai, bed_file ->
            def newMeta = meta.clone()
            newMeta.remove('bed_file')
            [newMeta, bam, bai, bed_file]
        }
        .set { ch_bams_bed_wgs }

    // Run mosdepth to get coverages
    MOSDEPTH(
        ch_bams_bed_targeted.mix(ch_bams_bed_wgs),
        fasta,
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.map { _meta, file -> file }.collect())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.map { _meta, file -> file }.collect())

    // Remove laneId, read_group, flowcellId, bed_file from the metadata to enable sample based grouping
    FASTP.out.json
        .map { meta, json ->
            def newMeta = meta.clone()
            newMeta.remove('laneId')
            newMeta.remove('read_group')
            newMeta.remove('flowcellId')
            newMeta.remove('bed_file')
            newMeta.remove('runId')
            [newMeta + [id: newMeta.sample], json]
        }
        .set { ch_fastp_mosdepth }

    // Remove bed_file from the metadata to enable sample based grouping - this result is coming from alignments in the samplesheet
    FASTQ_ALIGN_BWA_MARKDUPLICATES.out.jsonstats
        .map { meta, json ->
            def newMeta = meta.clone()
            newMeta.remove('bed_file')
            newMeta.remove('runId')
            [newMeta + [id: newMeta.sample], json]
        }
        .set { ch_fastp_mosdepth_aligned }

    // Collect the results for comparison
    MOSDEPTH.out.summary_txt
        .join(MOSDEPTH.out.regions_bed, by: 0)
        .join(ch_fastp_mosdepth.mix(ch_fastp_mosdepth_aligned).groupTuple(), by: 0)
        .set { ch_fastp_mosdepth_merged }

    // Compare coverage with thresholds: writing the results file
    // input: FASTP Q30 ratio + mosdepth all genes + mosdepth target genes
    COMPARE_THRESHOLD(
        ch_fastp_mosdepth_merged
    )
    ch_versions = ch_versions.mix(COMPARE_THRESHOLD.out.versions)


    // Merge compare_threshold results for a final report
    MERGE_REPORTS(
        COMPARE_THRESHOLD.out.result_csv.collect()
    )

    ch_multiqc_files = ch_multiqc_files.mix(MERGE_REPORTS.out.multiqc)
    ch_versions = ch_versions.mix(MERGE_REPORTS.out.versions)

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'pipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
