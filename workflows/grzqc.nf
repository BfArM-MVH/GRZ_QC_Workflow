/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_grzqc_pipeline'
include { CAT_FASTQ                       } from '../modules/nf-core/cat/fastq'
include { FASTQC                          } from '../modules/nf-core/fastqc'
include { FASTP                           } from '../modules/nf-core/fastp'
include { MULTIQC                         } from '../modules/nf-core/multiqc'
include { CONVERT_BED_CHROM               } from '../modules/local/convert_bed_chrom'
include { COMPARE_THRESHOLD               } from '../modules/local/compare_threshold'
include { MERGE_REPORTS                   } from '../modules/local/merge_reports'
include { BWAMEM2_INDEX                   } from '../modules/nf-core/bwamem2/index'
include { MOSDEPTH                        } from '../modules/nf-core/mosdepth'
include { SAMTOOLS_FAIDX                  } from '../modules/nf-core/samtools/faidx'
include { SAVE_REFERENCE                  } from '../modules/local/save_reference'
include { FASTQ_ALIGN_BWA_MARKDUPLICATES  } from '../subworkflows/local/fastq_align_bwa_markduplicates'

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

    // match fa and fasta exntensions
    def fastaExts = [ 'genome.fa', 'genome.fasta' ]
    def faiExts   = [ 'genome.fa.fai', 'genome.fasta.fai' ]

    // create reference channels
    if( params.reference_path ) {
        
        fasta = ch_genome
            .map { genome ->
            def candidates = fastaExts.collect { ext -> file("${params.reference_path}/${genome}/${ext}")}
            def f = candidates.find { it.exists() }            
            if( !f ) 
                error "Reference FASTA missing: $f"
            tuple( [ id: f.baseName ], f)
            }

        fai = ch_genome
            .map { genome ->
            def candidates = faiExts.collect { ext -> file("${params.reference_path}/${genome}/${ext}")}
            def f = candidates.find { it.exists() }            
            if( !f ) 
                error "Reference FASTA missing: $f"
            tuple( [ id: f.baseName ], f)
            }

        bwa = ch_genome
            .map { genome ->
            def f = file("${params.reference_path}/${genome}/bwamem2")
            if( !f.exists() ) error "BWA binary missing: $f"
            tuple('bwa', f)
            }

    }else{

        fasta = params.fasta ?
            // 1) user has supplied --fasta: load that single file
            Channel
            .fromPath(params.fasta, checkIfExists: true)
            .map { file -> tuple([ id: file.baseName ], file) }
        :
            // 2) no --fasta: pick default by genome
            ch_genome
            .flatMap { genome ->
                def defaultFasta = genome == 'GRCh38'
                    ? "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
                    : "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
                def f = file(defaultFasta)
                if( !f.exists() )    error "Default genome on ignomes s3 missing: $f"
                return f        
            }
            .map { file -> tuple([ id: file.baseName ], file) }
                    
        bwa         = params.bwa        ? Channel.fromPath(params.bwa ).map{ it -> [ [id:'bwa'], it ] }
                                        : Channel.empty()

        fai         = params.fai        ? Channel.fromPath(params.fai, checkIfExists: true).map{ file -> tuple([id: file.getSimpleName()], file) }
                                        : Channel.empty()
    }  

    // TARGET BED channel
    if( params.target ) {
        target = Channel.fromPath(params.target, checkIfExists: true)
    } else {
        target = ch_genome
                        .flatMap { genome ->
                        def defaultTargetCh = genome == 'GRCh38'
                            ? "${projectDir}/assets/default_files/hg38_440_omim_genes.bed"
                            : "${projectDir}/assets/default_files/hg19_439_omim_genes.bed"
                        def f = file(defaultTargetCh)
                        if( !f.exists() )    error "Default target BED missing: $f"
                        return f            
                        }
    }

    // CHROMOSOME‐NAME MAPPING channel
    if( params.mapping_chrom ) {
        mapping_chrom = Channel.fromPath(params.mapping_chrom, checkIfExists: true)
    } else {
        mapping_chrom = ch_genome
                            .flatMap { genome ->
                            def defaultMappingCh = genome == 'GRCh38'
                                ? "${projectDir}/assets/default_files/hg38_NCBI2UCSC.txt"
                                : "${projectDir}/assets/default_files/hg19_NCBI2UCSC.txt"
                            def f = file(defaultMappingCh)
                            if( !f.exists() )    error "Default mapping-chrom file missing: $f"
                            return f                                   
                            }
    }

    // create thresholds channels
    thresholds  = params.thresholds ? Channel.fromPath(params.thresholds, checkIfExists: true)
                                    : Channel.fromPath("${projectDir}/assets/default_files/thresholds.json")

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

    ch_bams = Channel.empty()

    if (!params.bwa && !params.reference_path){

        BWAMEM2_INDEX(
            fasta)

        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions.first())
        bwa = BWAMEM2_INDEX.out.index
    }

    if (!params.fai && !params.reference_path)
    {
        // Create sequence fai files for picard markduplicates
        //
        // MODULE: SAMTOOLS_FAIDX
        //
        SAMTOOLS_FAIDX(
            fasta,
            [[],[]],
            false)

        fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
    }

    if (params.save_reference){
        //save reference for the first run
        //
        // MODULE: SAVE_REFERENCE
        //
        SAVE_REFERENCE(
            fasta,
            fai,
            bwa,
            ch_genome
        )
    }

    //
    // SUBWORKFLOW: FASTQ_ALIGN_BWA_MARKDUPLICATES
    //

    // alignment analysis and markduplicates
    FASTQ_ALIGN_BWA_MARKDUPLICATES (
        ch_cat_fastq,
        bwa,
        true,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.versions.first())

    ch_bams =  FASTQ_ALIGN_BWA_MARKDUPLICATES.out.bam.join(FASTQ_ALIGN_BWA_MARKDUPLICATES.out.bai, by:0)

    // prepare mosdepth inputs
    // Prepare bed files
    // for WGS: defined bed files of ~400 genes
    // for WES and panel: extract bed from samplesheet and
    //      do the conversion with the correct mapping file (if it is NCBI format, covert it to UCSC).
    // for both cases different file required for hg38 and hg19
    ch_samplesheet_target = ch_samplesheet.filter { meta, fastqs, bed_file ->
        meta.libraryType in ["wes", "panel", "wes_lr", "panel_lr"]
    }
    ch_samplesheet_wgs = ch_samplesheet.filter { meta, fastqs, bed_file ->
        meta.libraryType in ["wgs", "wgs_lr"]
    }

    // --- For target samples: prepare bed for conversion ---
    ch_samplesheet_target.map { meta, fastqs, bed_file ->
        def bed = bed_file.first()
        return tuple(meta, bed)
    }.set{ch_bed_for_conversion_target}

    //
    // MODULE: CONVERT_BED_CHROM
    //
    // for WES and panel, run the conversion process: if the bed file has NCBI-style names, they will be converted.
    CONVERT_BED_CHROM (
        ch_bed_for_conversion_target,
        mapping_chrom
    )
    ch_converted_bed = CONVERT_BED_CHROM.out.converted_bed
    ch_versions = ch_versions.mix(CONVERT_BED_CHROM.out.versions)

    // Prepare mosdepth inputs for target samples
    ch_mosdepth_input_target = ch_samplesheet_target
        .map { meta, fastqs, bed_file -> [meta] }
        .join(ch_bams, by: 0)
        .join(ch_converted_bed, by: 0)


    // --- Prepare mosdepth inputs for WGS samples ---
    ch_mosdepth_input_wgs = ch_samplesheet_wgs
        .map { meta, fastqs, bed_file -> [meta] }
        .join(ch_bams, by: 0)
        .combine(target)

    // --- Merge the two sets of mosdepth inputs ---
    ch_mosdepth_input_target
        .mix(ch_mosdepth_input_wgs)
        .set{ch_mosdepth_input}

    //
    // MODULE: MOSDEPTH
    //
    MOSDEPTH(
        ch_mosdepth_input,
        fasta,
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.map{meta, file -> file}.collect())

    // collect the results for comparison
    MOSDEPTH.out.summary_txt.join(MOSDEPTH.out.regions_bed, by:0)
        .join(FASTP.out.json, by:0)
        .map{meta, summary, json, bed -> tuple(meta, json, summary,bed)}
        .set{ch_fastp_mosdepth_merged}

    //
    // MODULE:COMPARE_THRESHOLD
    //
    //Compare coverage with thresholds: writing the results file
    // input: FASTP Q30 ratio + mosdepth all genes + mosdepth target genes
    //
    COMPARE_THRESHOLD(
        ch_fastp_mosdepth_merged,
        thresholds
    )
    ch_versions = ch_versions.mix(COMPARE_THRESHOLD.out.versions)

    //
    // MODULE: MERGE_REPORTS
    //
    MERGE_REPORTS(
        COMPARE_THRESHOLD.out.result_csv.collect())

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
