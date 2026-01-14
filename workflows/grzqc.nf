/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                              } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQC                                        } from '../modules/nf-core/fastqc'
include { FASTP                                         } from '../modules/nf-core/fastp'
include { FASTPLONG                                     } from '../modules/local/fastplong'
include { MULTIQC                                       } from '../modules/nf-core/multiqc'
include { CONVERT_BED_CHROM                             } from '../modules/local/convert_bed_chrom'
include { COMPARE_THRESHOLD                             } from '../modules/local/compare_threshold'
include { MOSDEPTH                                      } from '../modules/nf-core/mosdepth'
include { FASTQ_ALIGN_BWA_MARKDUPLICATES                } from '../subworkflows/local/fastq_align_bwa_markduplicates'
include { SAMTOOLS_FASTQ                                } from '../modules/nf-core/samtools/fastq/main'
include { ALIGN_MERGE_LONG                              } from '../subworkflows/local/align_merge_long'
include { PREPARE_REFERENCES                            } from '../subworkflows/local/prepare_references'
include { MERGE_REPORTS as MERGE_REPORTS_DEDUPLICATED   } from '../modules/local/merge_reports'
include { MERGE_REPORTS as MERGE_REPORTS_UNDEDUPLICATED } from '../modules/local/merge_reports'

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

    // split rows into those starting from reads and those starting from alignments
    ch_samplesheet
        .branch { meta, reads, alignment ->
            reads: reads.size() > 0
            return [meta, reads]
            alignments: alignment.size() > 0
            return [meta, alignment]
        }
        .set { samplesheet_ch }

    // split rows starting from reads into short and long reads
    samplesheet_ch.reads
        .branch { meta, _reads ->
            lng: meta.libraryType.endsWith('_lr')
            srt: true
        }
        .set { samplesheet_ch_reads }

    PREPARE_REFERENCES(
        samplesheet_ch_reads.srt,
        samplesheet_ch_reads.lng,
        ch_genome,
    )
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    def fasta = PREPARE_REFERENCES.out.fasta
    def fai = PREPARE_REFERENCES.out.fai
    def bwa = PREPARE_REFERENCES.out.bwa
    def mmi = PREPARE_REFERENCES.out.mmi

    // split rows starting from alignments into short and long read alignments
    samplesheet_ch.alignments
        .branch { meta, _alignment ->
            lng: meta.libraryType.endsWith('_lr')
            srt: true
        }
        .set { samplesheet_ch_alignments }

    // split long read rows depending on if they are BAM (PacBio) or FASTQ (usually Nanopore)
    samplesheet_ch_reads.lng
        .branch { _meta, reads ->
            bam: reads.first().getExtension() == "bam"
            fastq: reads.first().getExtension() == "gz"
        }
        .set { samplesheet_ch_reads_lng }

    // Convert read BAMs to FASTQs
    def samtools_fastq_interleave = false
    SAMTOOLS_FASTQ(
        samplesheet_ch_reads_lng.bam,
        samtools_fastq_interleave,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

    // Run FASTQC on short + long FASTQ files - per lane
    // using the 'other' output because long-reads won't be paired-end
    FASTQC(
        samplesheet_ch_reads.srt.mix(samplesheet_ch_reads_lng.fastq).mix(SAMTOOLS_FASTQ.out.other)
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // Run FASP on FASTQ files - per lane
    save_trimmed_fail = false
    save_merged = false
    FASTP(
        samplesheet_ch_reads.srt.map{ meta, reads -> [meta, reads, []] },
        false,
        save_trimmed_fail,
        save_merged
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { _meta, json -> json })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect { _meta, html -> html })
    ch_versions = ch_versions.mix(FASTP.out.versions)

    FASTPLONG(
        samplesheet_ch_reads_lng.fastq.mix(SAMTOOLS_FASTQ.out.other),
        [],
        false,
        false,
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG.out.json.collect { _meta, json -> json })
    ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG.out.html.collect { _meta, html -> html })
    ch_versions = ch_versions.mix(FASTPLONG.out.versions)

    // align FASTQs per lane, merge, and sort
    FASTQ_ALIGN_BWA_MARKDUPLICATES(
        FASTP.out.reads,
        samplesheet_ch_alignments.srt,
        bwa,
        true,
        fasta,
        fai,
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.stat.collect { _meta, file -> file })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.flagstat.collect { _meta, file -> file })

    ch_reads_long_trimmed = FASTPLONG.out.reads.map { meta, reads ->
        def sequencer_manufacturer = meta.sequencer.toLowerCase()
        def is_pacbio = sequencer_manufacturer.contains("pacbio") || sequencer_manufacturer.contains("pacific biosciences")
        def mm2_preset = is_pacbio ? "map-hifi" : "map-ont"

        [meta + [mm2_preset: mm2_preset], reads]
    }
    ALIGN_MERGE_LONG(
        ch_reads_long_trimmed,
        samplesheet_ch_alignments.lng,
        mmi,
        fasta,
        fai,
    )
    ch_versions = ch_versions.mix(ALIGN_MERGE_LONG.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_MERGE_LONG.out.stat.collect { _meta, file -> file })
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_MERGE_LONG.out.flagstat.collect { _meta, file -> file })

    ch_bams = FASTQ_ALIGN_BWA_MARKDUPLICATES.out.bam.join(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.bai, by: 0).mix(ALIGN_MERGE_LONG.out.bam.join(ALIGN_MERGE_LONG.out.bai, by: 0))

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
        .mix(FASTPLONG.out.json)
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
            [newMeta + [id: newMeta.sample], json]
        }
        .set { ch_fastp_mosdepth_aligned }

    // duplicate fastp outputs for dedup and non-dedup comparison

    ch_fastp_fordedup = ch_fastp_mosdepth.mix(ch_fastp_mosdepth_aligned)
    ch_fastp_skipdedup = ch_fastp_mosdepth.mix(ch_fastp_mosdepth_aligned)

    ch_fastp_fordedup.map { meta, json -> [meta + [is_deduplicated: true], json] }.set { ch_fastp_mosdepth_duplicated }

    ch_fastp_skipdedup.map { meta, json -> [meta + [is_deduplicated: false], json] }.set { ch_fastp_mosdepth_nodedup }

    // Collect the results for comparison
    MOSDEPTH.out.summary_txt
        .join(MOSDEPTH.out.regions_bed, by: 0)
        .join(ch_fastp_mosdepth_duplicated.mix(ch_fastp_mosdepth_nodedup).groupTuple(), by: 0)
        .set { ch_fastp_mosdepth_merged }

    // Compare coverage with thresholds: writing the results file
    // input: FASTP Q30 ratio + mosdepth all genes + mosdepth target genes

    COMPARE_THRESHOLD(
        ch_fastp_mosdepth_merged
    )
    ch_versions = ch_versions.mix(COMPARE_THRESHOLD.out.versions)

    COMPARE_THRESHOLD.out.result_csv
        .branch { meta, _file ->
            deduplicated: meta.is_deduplicated == true
            undeduplicated: meta.is_deduplicated == false
        }
        .set { ch_compare_threshold_results }

    // Merge compare_threshold results for a final report
    MERGE_REPORTS_DEDUPLICATED(
        ch_compare_threshold_results.deduplicated.collect { _meta, file -> file }
    )
    ch_multiqc_files = ch_multiqc_files.mix(MERGE_REPORTS_DEDUPLICATED.out.multiqc)
    ch_versions = ch_versions.mix(MERGE_REPORTS_DEDUPLICATED.out.versions)

    MERGE_REPORTS_UNDEDUPLICATED(
        ch_compare_threshold_results.undeduplicated.collect { _meta, file -> file }
    )

    ch_multiqc_files = ch_multiqc_files.mix(MERGE_REPORTS_UNDEDUPLICATED.out.multiqc)
    ch_versions = ch_versions.mix(MERGE_REPORTS_UNDEDUPLICATED.out.versions)

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
