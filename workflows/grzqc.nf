/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_grzqc_pipeline'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq'
include { FASTQC                 } from '../modules/nf-core/fastqc'
include { FASTP                  } from '../modules/nf-core/fastp'
include { MULTIQC                } from '../modules/nf-core/multiqc'
include { COMPARE_THRESHOLD      } from '../modules/local/compare_threshold'
include { MERGE_REPORTS          } from '../modules/local/merge_reports'
include { FASTQ_ALIGN_BWA as FASTQ_ALIGN_BWA_HG38 } from '../subworkflows/nf-core/fastq_align_bwa'
include { FASTQ_ALIGN_BWA as FASTQ_ALIGN_BWA_HG37 } from '../subworkflows/nf-core/fastq_align_bwa'
include { BWAMEM2_INDEX  as  BWAMEM2_INDEX_HG38   } from '../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_INDEX  as  BWAMEM2_INDEX_HG37   } from '../modules/nf-core/bwamem2/index/main'
include { MOSDEPTH_MOSDEPTH_TARGET as MOSDEPTH_MOSDEPTH_TARGET_HG38 } from '../subworkflows/local/mosdepth_mosdepth_target'
include { MOSDEPTH_MOSDEPTH_TARGET as MOSDEPTH_MOSDEPTH_TARGET_HG37 } from '../subworkflows/local/mosdepth_mosdepth_target'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GRZQC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // create reference channels
    fasta_37       = Channel.fromPath(params.fasta_37, checkIfExists: true)
                    .map{ file -> tuple([id: file.getSimpleName()], file) }.collect()
    fasta_38        = Channel.fromPath(params.fasta_38, checkIfExists: true)
                    .map{ file -> tuple([id: file.getSimpleName()], file) }.collect()

    rep_genes_38  = Channel.fromPath(params.rep_genes_hg38, checkIfExists: true).collect()

    rep_genes_37  = Channel.fromPath(params.rep_genes_hg19, checkIfExists: true).collect()

    bwa_index_37    = params.bwa_index_37   ? Channel.fromPath(params.bwa_index_37, checkIfExists: true).map{ file -> tuple([id: file.getSimpleName()], file) }.collect()
                                            : Channel.empty()

    bwa_index_38    = params.bwa_index_38   ? Channel.fromPath(params.bwa_index_38, checkIfExists: true).map{ file -> tuple([id: file.getSimpleName()], file) }.collect()
                                            : Channel.empty()

    thresholds  = Channel.fromPath("${projectDir}/assets/thresholds.json", checkIfExists: true).collect()

    ch_samplesheet.view()
    // Creating merge fastq channel from samplesheet
    ch_samplesheet.map{
        meta, fastqs, bed_file -> tuple(meta, fastqs)
    }.branch {
            meta, fastqs ->
                single_sample_paired_end  : fastqs.size() < 2 | fastqs.size() == 2
                    return [ meta, fastqs.flatten() ]
                multiple_sample_paired_end: fastqs.size() > 2
                    return [ meta, fastqs.flatten() ]
        }.set{
            ch_fastqs
        }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastqs.multiple_sample_paired_end
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    CAT_FASTQ.out.reads
        .mix(ch_fastqs.single_sample_paired_end)
        .set { ch_cat_fastq }


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_cat_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: FASTP
    //
    save_trimmed_fail = false
    save_merged = false
    FASTP(
        ch_cat_fastq,
        [], // we are not using any adapter fastas at the moment
        false, // we don't use discard_trimmed_pass at the moment
        save_trimmed_fail,
        save_merged
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{ meta, json -> json })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect{ meta, html -> html })
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    // branch out samples per reference
    ch_cat_fastq.branch{
        def meta = it[0]
        hg38: meta.reference == "GRCh38" || meta.reference == "hg38"
        hg37: meta.reference == "GRCh37" || meta.reference == "hg19"
        other:false
    }.set{input_samples}

    ch_bams = Channel.empty()
    // bwa index creation might be dropped
    if (!params.bwa_index_38){
        BWAMEM2_INDEX_HG38(fasta_38)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_HG38.out.versions.first())
        bwa_index_38 = BWAMEM2_INDEX_HG38.out.index
    }

    // hg38 alignment analysis
    FASTQ_ALIGN_BWA_HG38 (
        input_samples.hg38,
        bwa_index_38,
        true,
        fasta_38
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA_HG38.out.versions.first())
    ch_bams = ch_bams
                .mix(FASTQ_ALIGN_BWA_HG38.out.bam
                    .join(FASTQ_ALIGN_BWA_HG38.out.bai, by:0))

    // bwa index creation might be dropped
    if (!params.bwa_index_37){
        BWAMEM2_INDEX_HG37(fasta_37)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX_HG38.out.versions.first())
        bwa_index_37 = BWAMEM2_INDEX_HG37.out.index
    }
    // hg37 alignment analysis
    FASTQ_ALIGN_BWA_HG37 (
        input_samples.hg37,
        bwa_index_37,
        true,
        fasta_37
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA_HG37.out.versions.first())
    ch_bams = ch_bams
                .mix(FASTQ_ALIGN_BWA_HG37.out.bam
                    .join(FASTQ_ALIGN_BWA_HG37.out.bai, by:0))

    // prepare mosdepth inputs
    ch_samplesheet
        .map{meta, fastqs, bed_file -> tuple(meta, bed_file.first())}
        .join(ch_bams, by:0)
        .map{meta, bed_file, bam, bai -> tuple(meta, bam, bai, bed_file )}
        .branch{
            def meta = it[0]
            hg38: meta.reference == "GRCh38" || meta.reference == "hg38"
            hg37: meta.reference == "GRCh37" || meta.reference == "hg19"
            other:false
        }
        .set{ch_mosdepth_input}

    // hg38 mosdepth analysis
    MOSDEPTH_MOSDEPTH_TARGET_HG38(
        ch_mosdepth_input.hg38,
        fasta_38,
        rep_genes_38
    )
    ch_versions = ch_versions.mix(MOSDEPTH_MOSDEPTH_TARGET_HG38.out.ch_versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_MOSDEPTH_TARGET_HG38.out.ch_multiqc_files)


    // hg37 mosdepth analysis
    MOSDEPTH_MOSDEPTH_TARGET_HG37(
        ch_mosdepth_input.hg37,
        fasta_37,
        rep_genes_37
    )
    ch_versions = ch_versions.mix(MOSDEPTH_MOSDEPTH_TARGET_HG37.out.ch_versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_MOSDEPTH_TARGET_HG37.out.ch_multiqc_files)


    // collect the results for comparison

    MOSDEPTH_MOSDEPTH_TARGET_HG37.out.summary_txt.mix(
        MOSDEPTH_MOSDEPTH_TARGET_HG38.out.summary_txt
    ).set{ch_mosdepth_summary}

    MOSDEPTH_MOSDEPTH_TARGET_HG37.out.regions_bed.mix(
        MOSDEPTH_MOSDEPTH_TARGET_HG38.out.regions_bed
    ).set{ch_mosdepth_target_regions}

    ch_mosdepth_summary.join(ch_mosdepth_target_regions, by:0)
        .join(FASTP.out.json, by:0)
        .map{meta, summary, json, bed -> tuple(meta, json, summary,bed)}
        .set{ch_fastp_mosdepth_merged}

    //
    // MODULE: Compare coverage with thresholds: writing the results file
    // input: FASTP Q30 ratio + mosdepth all genes + mosdepth target genes
    //
    COMPARE_THRESHOLD(
        ch_fastp_mosdepth_merged,
        thresholds
    )
    ch_versions = ch_versions.mix(COMPARE_THRESHOLD.out.versions)

    //
    // MODULE: Merge Reports
    //
    MERGE_REPORTS(COMPARE_THRESHOLD.out.result_csv.collect())
    ch_versions = ch_versions.mix(MERGE_REPORTS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
