//
// Alignment with BWA, merge with Samtools and markduplicates with sambamba
//

include { BWAMEM2_MEM              } from '../../../modules/nf-core/bwamem2/mem/main'
include { BAM_INDEX_STATS_SAMTOOLS } from '../../local/bam_index_stats_samtools/main'
include { FASTQ_SORT               } from '../../../modules/local/fastq_sort/main'
include { SAMTOOLS_MERGE           } from '../../../modules/nf-core/samtools/merge/main'
include { SAMBAMBA_MARKDUP         } from '../../../modules/nf-core/sambamba/markdup/main'

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
            tuple(meta, meta.fastp_json)
        }
        .set { jsonstats }

    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    BAM_INDEX_STATS_SAMTOOLS(
        SAMTOOLS_MERGE.out.bam.mix(ch_alignments_newMeta),
        ch_fasta,
    )

    ch_versions = ch_versions.mix(BAM_INDEX_STATS_SAMTOOLS.out.versions)
    ch_bam = BAM_INDEX_STATS_SAMTOOLS.out.bam
    ch_bai = BAM_INDEX_STATS_SAMTOOLS.out.bai

    if (!params.skip_markdup) {
        // Mark duplicates with sambamba
        SAMBAMBA_MARKDUP(
            BAM_INDEX_STATS_SAMTOOLS.out.bam
        )
        ch_versions = ch_versions.mix(SAMBAMBA_MARKDUP.out.versions)
        ch_bam = SAMBAMBA_MARKDUP.out.bam
        // channel: [ val(meta), path(bam) ]
        ch_bai = SAMBAMBA_MARKDUP.out.bai
    }

    emit:
    bam       = ch_bam // channel: [ val(meta), path(bam) ]
    bai       = ch_bai // channel: [ val(meta), path(bai) ]
    flagstat  = BAM_INDEX_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    stat      = BAM_INDEX_STATS_SAMTOOLS.out.stats // channel: [ val(meta), path(stats) ]
    jsonstats // channel: [ val(meta), path(json) ]
    versions  = ch_versions // channel: [ path(versions.yml) ]
}
