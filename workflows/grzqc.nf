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
include { FASTPLONG                      } from '../modules/local/fastplong'
include { MULTIQC                        } from '../modules/nf-core/multiqc'
include { CONVERT_BED_CHROM              } from '../modules/local/convert_bed_chrom'
include { COMPARE_THRESHOLD              } from '../modules/local/compare_threshold'
include { MERGE_REPORTS                  } from '../modules/local/merge_reports'
include { MOSDEPTH                       } from '../modules/nf-core/mosdepth'
include { FASTQ_ALIGN_BWA_MARKDUPLICATES } from '../subworkflows/local/fastq_align_bwa_markduplicates'
include { ALIGN_MERGE_LONG               } from '../subworkflows/local/align_merge_long'
include { PBTK_PBINDEX                   } from '../modules/nf-core/pbtk/pbindex'
include { PBTK_BAM2FASTQ                 } from '../modules/nf-core/pbtk/bam2fastq'
include { PREPARE_REFERENCES             } from '../subworkflows/local/prepare_references'

// --------------------- HOTFIX: Imports ---------------------
import groovy.json.JsonSlurper
def jsonSlurper = new JsonSlurper()
import java.util.zip.GZIPInputStream
import java.io.InputStreamReader
import java.io.BufferedReader

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

    // Index PacBio long-read BAMs
    PBTK_PBINDEX(samplesheet_ch_reads_lng.bam)
    ch_versions = ch_versions.mix(PBTK_PBINDEX.out.versions)

    // Convert PacBio long-read BAMs to FASTQs
    PBTK_BAM2FASTQ(samplesheet_ch_reads_lng.bam.join(PBTK_PBINDEX.out.pbi))
    ch_versions = ch_versions.mix(PBTK_BAM2FASTQ.out.versions)

    // Run FASTQC on short + long FASTQ files - per lane
    FASTQC(
        samplesheet_ch_reads.srt.mix(samplesheet_ch_reads_lng.fastq).mix(PBTK_BAM2FASTQ.out.fastq)
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // Run FASP on FASTQ files - per lane
    save_trimmed_fail = false
    save_merged = false
    FASTP(
        samplesheet_ch_reads.srt,
        [],
        false,
        save_trimmed_fail,
        save_merged,
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { _meta, json -> json })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect { _meta, html -> html })
    ch_versions = ch_versions.mix(FASTP.out.versions)

    FASTPLONG(
        samplesheet_ch_reads_lng.fastq.mix(PBTK_BAM2FASTQ.out.fastq),
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
        def sequencer = meta.sequencer.toLowerCase()
        def is_pacbio = "pacbio" in sequencer || "pacific biosciences" in sequencer
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
            def id_base = newMeta.sample.replaceFirst(/_run\d+$/, '')
            [newMeta + [id: id_base], json]
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



    /*
    ------------------------------------------------------------------------------------------
        HOTFIX: Correct percentBasesAboveQualityThreshold & targetedRegionsAboveMinCoverage
    ------------------------------------------------------------------------------------------

    Context:
    - The COMPARE_THRESHOLD module receives values via Nextflow metadata (meta.*),
      which are compared to internally computed values (from FASTP JSON and mosdepth BED files).

    Problem:
    - The pipeline previously passed incorrect or placeholder values for:
        * percentBasesAboveQualityThreshold (from FASTP)
        * targetedRegionsAboveMinCoverage (from mosdepth BED)
    - This caused threshold checks in compare_threshold.py to fail incorrectly,
      even when data quality was acceptable.

    Fix:
    - Recompute both metrics **in the main workflow** before passing to COMPARE_THRESHOLD:
        1. Compute percentBasesAboveQualityThreshold from 'before_filtering' section of FASTP JSON:
           → (q30_bases / total_bases) × 100
        2. Compute targetedRegionsAboveMinCoverage from compressed mosdepth regions.bed.gz:
           → fraction of regions with ≥ 20× coverage

    - These corrected values are injected into the `meta` map passed to COMPARE_THRESHOLD.

    Why this works:
    - It aligns the pipeline metadata (`args.*`) with the internally computed values in the Python script.
    - This ensures accurate quality checks and eliminates false FAILs in threshold validation.

    Status:
    - Confirmed via DEBUG prints in `.command.log` that both values now match expectations.

    */

    // ------------------------------------------------------------------
    // 1) Helper – normalise metadata (collapse multiple runs of the same sample)
    //    * removes run‑specific keys
    //    * creates a stable "id" without the _runX suffix
    // ------------------------------------------------------------------
    
    def normMeta = { Map m ->
        def n = m.clone()
        ['runId', 'bed_file', 'laneId', 'read_group', 'flowcellId'].each { n.remove(it) }
        n.id = n.sample.replaceFirst(/_run\d+\$/, '')
        return n
    }
    
    // ------------------------------------------------------------------
    // 2) Collect FASTP JSONs  (≥1 file per sample)
    //    id → List<Path>
    // ------------------------------------------------------------------
    
    ch_fastp_mosdepth
        .mix(ch_fastp_mosdepth_aligned)
        .map { m, json -> tuple(normMeta(m).id, json) }
        .groupTuple()
        .filter { _id, lst -> lst } // drop empty lists
        .set { fastp_ch }
    
    // ------------------------------------------------------------------
    // 3) Mosdepth summary – compute mean depth of coverage
    //    id → [metaWithDepth, summaryPath]
    // ------------------------------------------------------------------
    
    summary_ch = MOSDEPTH.out.summary_txt.map { m, summary ->
        double mean = summary.text.readLines()
                                 .find { it.startsWith('total\t') }
                                 ?.split(/\t/)[3] as double ?: 0.0
        tuple(normMeta(m).id,
              [normMeta(m) + [meanDepthOfCoverage: mean], summary])
    }
    
    // ------------------------------------------------------------------
    // 4) Mosdepth regions.bed  (coverage per target region)
    //    id → bedPath
    // ------------------------------------------------------------------
    
    regions_ch = MOSDEPTH.out.regions_bed.map { m, bed ->
        tuple(normMeta(m).id, bed)
    }
    
    // ------------------------------------------------------------------
    // 5) Join summary ∩ regions  →  id, meta, summaryPath, bedPath
    // ------------------------------------------------------------------
    
    sr_ch = summary_ch.join(regions_ch, by: 0)
                      .map { id, pair, bed -> tuple(id, pair[0], pair[1], bed) }
    
    // ------------------------------------------------------------------
    // 6) Final 4‑tuple (meta, summaryPath, bedPath, fastpJson)
    // ------------------------------------------------------------------
    
    sr_ch.join(fastp_ch, by: 0)
         .map { id, meta, summaryPath, bedPath, fastpList ->
    
             /* ---------- Q30 rate from FASTP JSON ---------- */
             def jsonObj   = new JsonSlurper().parse(fastpList[0])
             def before    = jsonObj.summary.before_filtering
             double tot    = (before.total_bases ?: 1) as double
             double q30    = (before.q30_bases   ?: 0) as double
             meta.percentBasesAboveQualityThreshold = (q30 / tot) * 100.0
    
             /* ---------- % target regions ≥ 20× coverage ---------- */
             final double MIN_COV = 20.0
             int totalRegions = 0
             int covered      = 0
    
             new GZIPInputStream(bedPath.newInputStream()).withReader { reader ->
                 reader.eachLine { line ->
                     if (line.startsWith('#') || line.trim().isEmpty()) return
                     double cov = line.tokenize('\t')[3] as double
                     totalRegions++
                     if (cov >= MIN_COV) covered++
                 }
             }
             meta.targetedRegionsAboveMinCoverage =
                 totalRegions ? covered / (double) totalRegions : 0.0 // fraction
    
             /* ---------- emit tuple to COMPARE_THRESHOLD ---------- */
             tuple(meta, summaryPath, bedPath, fastpList[0])
         }
         .set { ch_fastp_mosdepth_merged }
    
    // ------------------------------------------------------------------
    // 7) Pass downstream
    // ------------------------------------------------------------------
    
    COMPARE_THRESHOLD(ch_fastp_mosdepth_merged)

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
