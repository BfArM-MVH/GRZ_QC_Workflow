//
// Align with minimap2 and merge with samtools
//

include { MINIMAP2_ALIGN           } from '../../../modules/nf-core/minimap2/align/main'
include { BAM_INDEX_STATS_SAMTOOLS } from '../../local/bam_index_stats_samtools/main'
include { SAMTOOLS_MERGE           } from '../../../modules/nf-core/samtools/merge/main'
include { CALCULATE_BASEQUALITY    } from '../../../modules/local/calculate_basequality/main'

workflow ALIGN_MERGE_LONG {
    take:
    ch_reads      // channel (optional): [ val(meta), [ path(reads) ] ]
    ch_alignments // channel (optional): [ val(meta), path(alignments) ]
    ch_index      // channel (mandatory): [ val(meta2), path(index) ]
    ch_fasta      // channel (mandatory) : [ val(meta3), path(fasta) ]
    ch_fai        // channel (mandatory) : [ val(meta4), path(fai) ]

    main:
    ch_versions = Channel.empty()

    // Map reads with minimap2 - per lane
    // (nf-core module always sorts BAM output)
    def mm2_output_bam = true
    def mm2_custom_bam_index_ext = ''
    def mm2_use_paf_cigar_fmt = false
    // work around BAM limit to CIGAR length which can be exceeded with long reads
    // https://github.com/lh3/minimap2#working-with-65535-cigar-operations
    def mm2_write_long_cigar_to_bam_tag = true
    MINIMAP2_ALIGN(
        ch_reads,
        ch_index,
        mm2_output_bam,
        mm2_custom_bam_index_ext,
        mm2_use_paf_cigar_fmt,
        mm2_write_long_cigar_to_bam_tag,
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    // Remove laneId, read_group, flowcellId from the metadata to enable sample based grouping
    ch_bams = MINIMAP2_ALIGN.out.bam
        .map { meta, bam ->
            def newMeta = meta.clone()
            newMeta.remove('laneId')
            newMeta.remove('read_group')
            newMeta.remove('flowcellId')
            newMeta.remove('runId')
            newMeta.remove('mm2_preset')
            [newMeta + [id: newMeta.sample], bam]
        }
        .groupTuple()
    // Merge alignments from different lanes
    SAMTOOLS_MERGE(
        ch_bams,
        ch_fasta,
        ch_fai,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)


    // run calculate_basequality.py on the alogned bam file (from samplesheet)
    CALCULATE_BASEQUALITY(
        ch_alignments
    )
    ch_versions = ch_versions.mix(CALCULATE_BASEQUALITY.out.versions)

    // Merge aligned bams with the alignments coming from samplesheet

    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    BAM_INDEX_STATS_SAMTOOLS(
        SAMTOOLS_MERGE.out.bam.mix(ch_alignments),
        ch_fasta,
    )

    ch_versions = ch_versions.mix(BAM_INDEX_STATS_SAMTOOLS.out.versions)

    emit:
    bam       = BAM_INDEX_STATS_SAMTOOLS.out.bam // channel: [ val(meta), path(bam) ]
    bai       = BAM_INDEX_STATS_SAMTOOLS.out.bai // channel: [ val(meta), path(bai) ]
    flagstat  = BAM_INDEX_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    stat      = BAM_INDEX_STATS_SAMTOOLS.out.stats // channel: [ val(meta), path(stats) ]
    jsonstats = CALCULATE_BASEQUALITY.out.json // channel: [ val(meta), path(json) ]
    versions  = ch_versions // channel: [ path(versions.yml) ]
}
