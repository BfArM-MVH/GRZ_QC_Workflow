//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../nf-core/bam_stats_samtools/main'

workflow BAM_INDEX_STATS_SAMTOOLS {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    ch_fai // channel: [ val(meta), path(fai) ]

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_INDEX(
        ch_bam
    )

    ch_bam
        .join(SAMTOOLS_INDEX.out.index, by: [0], remainder: true)
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS(
        ch_bam_bai,
        ch_fasta.combine(ch_fai).map { meta, fasta, _meta2, fai -> [meta, fasta, fai] }.first(),
    )

    emit:
    bam      = ch_bam // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.index // channel: [ val(meta), [ bai/csi ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
