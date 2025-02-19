//
// MOSDEPTH analysis
//

include { MOSDEPTH                      } from '../../../modules/nf-core/mosdepth'
include { MOSDEPTH as MOSDEPTH_TARGET   } from '../../../modules/nf-core/mosdepth'

workflow MOSDEPTH_MOSDEPTH_TARGET {
    take:
    ch_mosdepth_input  // channel (mandatory): [ val(meta), bam, bai, bed ]
    ch_fasta           // channel (optional) : [ val(meta3), path(fasta) ]
    ch_rep_genes       // channel (optional) : [ val(meta4), path(bed) ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: MOSDEPTH, get average coverage
    //
    MOSDEPTH (
        ch_mosdepth_input,
        ch_fasta
    )
    summary_txt = MOSDEPTH.out.summary_txt
    ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.map{meta, file -> file}.collect())
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    //
    // MODULE: MOSDEPTH, get coverage on target genes
    //
    ch_mosdepth_input.map{
        meta, bam, bai, bed -> tuple(meta, bam, bai)
    }.combine(ch_rep_genes)
    .set { ch_mosdepth_target_genes_input }

    MOSDEPTH_TARGET  (
        ch_mosdepth_target_genes_input,
        ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_TARGET.out.summary_txt.map{meta, file -> file}.collect())
    regions_bed = MOSDEPTH_TARGET.out.regions_bed
    ch_versions = ch_versions.mix(MOSDEPTH_TARGET.out.versions.first())

    emit:
    summary_txt
    regions_bed
    ch_multiqc_files
    ch_versions                          // channel: [ path(versions.yml) ]
}
