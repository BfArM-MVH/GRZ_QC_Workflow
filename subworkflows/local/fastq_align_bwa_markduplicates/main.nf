//
// Alignment with BWA, merge with Samtools and markduplicates with sambamba
//

include { BWAMEM2_MEM              } from '../../../modules/nf-core/bwamem2/mem/main'
include { BAM_INDEX_STATS_SAMTOOLS } from '../../local/bam_index_stats_samtools/main'
include { FASTQ_SORT               } from '../../../modules/local/fastq_sort/main'
include { SAMTOOLS_MERGE           } from '../../../modules/nf-core/samtools/merge/main'
include { SAMBAMBA_MARKDUP         } from '../../../modules/nf-core/sambamba/markdup/main'
include { UMITOOLS_DEDUP           } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX           } from '../../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN_BWA_MARKDUPLICATES {
    take:
    ch_reads      // channel (optional): [ val(meta), [ path(reads) ] ]
    ch_alignments // channel (optional): [ val(meta), path(alignments) ]
    ch_index      // channel (mandatory): [ val(meta2), path(index) ]
    val_sort_bam  // boolean (mandatory): true or false
    ch_fasta      // channel (mandatory) : [ val(meta3), path(fasta) ]
    ch_fai        // channel (mandatory) : [ val(meta4), path(fai) ]

    main:
    ch_versions = Channel.empty()
    ch_bam = Channel.empty()
    ch_bai = Channel.empty()

    ch_alignment_reads = ch_reads

    if (params.sort_paired_fastq) {
        // Sort paired-end FASTQ files (required for bwa-mem)
        FASTQ_SORT(
            ch_reads
        )

        ch_versions = ch_versions.mix(FASTQ_SORT.out.versions)

        ch_alignment_reads = FASTQ_SORT.out.reads
    }

    // Map reads with BWA - per lane
    BWAMEM2_MEM(
        ch_alignment_reads,
        ch_index,
        ch_fasta,
        val_sort_bam,
    )

    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    // Remove laneId, read_group, flowcellId from the metadata to enable sample based grouping
    ch_bams = BWAMEM2_MEM.out.bam
        .map { meta, bam ->
            def newMeta = meta.clone()
            newMeta.remove('laneId')
            newMeta.remove('read_group')
            newMeta.remove('flowcellId')
            newMeta.remove('runId')
            [newMeta + [id: newMeta.sample], bam]
        }
        .groupTuple()
    // Merge alignments from different lanes
    SAMTOOLS_MERGE(
        ch_bams,
        ch_fasta,
        ch_fai,
        [[], []],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // make sure ch_alignment has the same metadata for downstream joins
    def ch_alignments_newMeta = ch_alignments.map { meta, bam ->
        def newMeta = meta.clone()
        newMeta.remove('runId')
        newMeta.remove('laneId')
        newMeta.remove('flowcellId')
        [newMeta + [id: newMeta.sample], bam]
    }

    ch_alignments_newMeta
        .map { meta, _bam ->
            def newMeta = meta.clone()
            newMeta.remove('umi_base_skip')
            newMeta.remove('umi_length')
            newMeta.remove('umi_location')
            newMeta.remove('umi_in_read_header')  
            tuple(newMeta, newMeta.fastp_json)      
        }
        .set { jsonstats }

    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    BAM_INDEX_STATS_SAMTOOLS(
        SAMTOOLS_MERGE.out.bam.mix(ch_alignments_newMeta),
        ch_fasta,
    )

    ch_versions = ch_versions.mix(BAM_INDEX_STATS_SAMTOOLS.out.versions)
    //ch_bam = ch_bam.mix(BAM_INDEX_STATS_SAMTOOLS.out.bam.map { meta, bam -> [meta + [is_deduplicated: false], bam] })
    //ch_bai = ch_bai.mix(BAM_INDEX_STATS_SAMTOOLS.out.bai.map { meta, bai -> [meta + [is_deduplicated: false], bai] })

    // Deduplicate UMI tagged BAM files with umitools dedup - not optional for UMIreads
    BAM_INDEX_STATS_SAMTOOLS.out.bam.join(BAM_INDEX_STATS_SAMTOOLS.out.bai)
        .branch { meta, bam, bai ->
            umi: (meta.umi_base_skip || meta.umi_length || meta.umi_location || meta.umi_in_read_header)
            standard: true
        }
        .set { ch_for_dedup }

    // remove umi metadata for umitools dedup to avoid confusion in output metadata
    def umi_bams = ch_for_dedup.umi.map { meta, bam, bai ->
        def newMeta = meta.clone()
        newMeta.remove('umi_base_skip')
        newMeta.remove('umi_length')
        newMeta.remove('umi_location')
        newMeta.remove('umi_in_read_header')
        [newMeta, bam, bai]
    }

    def standard_bams = ch_for_dedup.standard.map { meta, bam, bai ->
        def newMeta = meta.clone()
        newMeta.remove('umi_base_skip')
        newMeta.remove('umi_length')
        newMeta.remove('umi_location')
        newMeta.remove('umi_in_read_header')
        [newMeta, bam, bai]
    }

    // We report both duplicated and deduplicated results for no UMI tagged samples, to enable comparison of the effect of deduplication on the results. 
    // For UMI tagged samples we only report non-deduplicated results. 
    ch_bam = ch_bam.mix(standard_bams.map{ meta, bam, bai -> [meta + [is_deduplicated: false], bam] })
    ch_bai = ch_bai.mix(standard_bams.map { meta, bam, bai -> [meta + [is_deduplicated: false], bai] })

    UMITOOLS_DEDUP(
        umi_bams,
        true, // get_output_stats
    )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)

    SAMTOOLS_INDEX(
        UMITOOLS_DEDUP.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    ch_bam = ch_bam.mix(UMITOOLS_DEDUP.out.bam.map { meta, bam -> [meta + [is_deduplicated: true], bam] })
    ch_bai = ch_bai.mix(SAMTOOLS_INDEX.out.bai.map { meta, bai -> [meta + [is_deduplicated: true], bai] })

    if (!params.skip_markdup) {
        // Mark duplicates with sambamba
        SAMBAMBA_MARKDUP(
            standard_bams
        )
        ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)
        ch_bam = ch_bam.mix(SAMBAMBA_MARKDUP.out.bam.map { meta, bam -> [meta + [is_deduplicated: true], bam] })
        ch_bai = ch_bai.mix(SAMBAMBA_MARKDUP.out.bai.map { meta, bai -> [meta + [is_deduplicated: true], bai] })
    }

    emit:
    bam       = ch_bam // channel: [ val(meta), path(bam) ]
    bai       = ch_bai // channel: [ val(meta), path(bai) ]
    flagstat  = BAM_INDEX_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    stat      = BAM_INDEX_STATS_SAMTOOLS.out.stats // channel: [ val(meta), path(stats) ]
    jsonstats // channel: [ val(meta), path(json) ]
    umi_log   = UMITOOLS_DEDUP.out.log // channel: [ val(meta), path(log) ]
    versions  = ch_versions // channel: [ path(versions.yml) ]
}
