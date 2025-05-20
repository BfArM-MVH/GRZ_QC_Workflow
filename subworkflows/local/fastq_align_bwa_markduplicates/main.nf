//
// Alignment with BWA
//

include { BWAMEM2_MEM                   } from '../../../modules/nf-core/bwamem2/mem/main'
include { BAM_SORT_STATS_SAMTOOLS       } from '../../nf-core/bam_sort_stats_samtools/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../../../modules/nf-core/picard/addorreplacereadgroups/main'      
include { PICARD_MARKDUPLICATES         } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_MERGE                } from '../../../modules/nf-core/samtools/merge/main'

workflow FASTQ_ALIGN_BWA_MARKDUPLICATES {
    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index        // channel (mandatory): [ val(meta2), path(index) ]
    val_sort_bam    // boolean (mandatory): true or false
    ch_fasta        // channel (mandatory) : [ val(meta3), path(fasta) ]
    ch_fai          // channel (mandatory) : [ val(meta4), path(fai) ]

    main:
    ch_versions = Channel.empty()

    // Map reads with BWA - per lane
    BWAMEM2_MEM ( 
        ch_reads, 
        ch_index, 
        ch_fasta, 
        val_sort_bam )

    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    // Merge alignments from different lanes
    ch_bams = BWAMEM2_MEM.out.bam.map { meta, bam -> 
                                            def newMeta = meta.clone()
                                            newMeta.remove('laneId')
                                            newMeta.remove('read_group')
                                            [ newMeta, bam ]}.groupTuple()

    SAMTOOLS_MERGE(
        ch_bams,
        ch_fasta,
        ch_fai
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    BAM_SORT_STATS_SAMTOOLS (
        SAMTOOLS_MERGE.out.bam, 
        ch_fasta)

    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)


    // run picard markduplicates on aligned and merged files
    PICARD_MARKDUPLICATES (
        SAMTOOLS_MERGE.out.bam, 
        ch_fasta, 
        ch_fai)

    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    emit:
    bam      = PICARD_MARKDUPLICATES.out.bam        // channel: [ val(meta), path(bam) ]
    bai      = PICARD_MARKDUPLICATES.out.bai        // channel: [ val(meta), path(bai) ]
    cram     = PICARD_MARKDUPLICATES.out.cram       // channel: [ val(meta), path(cram) ]
    metrics  = PICARD_MARKDUPLICATES.out.metrics    // channel: [ val(meta), path(metrics) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    stat     = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    versions = ch_versions                          // channel: [ path(versions.yml) ]
}
