#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BfArM-MVH/GRZ_QC_Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/BfArM-MVH/GRZ_QC_Workflow
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GRZQC                   } from './workflows/grzqc'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_grzqc_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_grzqc_pipeline'
include { METADATA_TO_SAMPLESHEET } from './modules/local/metadata_to_samplesheet'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_GRZQC {

    take:
    samplesheet // channel: samplesheet created by METADATA_TO_SAMPLESHEET
    genome      // string: genome 

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GRZQC (
        samplesheet,
        genome
    )
    emit:
    multiqc_report = GRZQC.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // Check whether a submission_basepath or a samplesheet is provided
    //
    def submission_basepath = (params.submission_basepath && params.submission_basepath.trim() != 'null') ? params.submission_basepath : null
    def input_samplesheet = (params.input && params.input.trim() != 'null') ? params.input : null

    if( submission_basepath && input_samplesheet ) {
        error "Please provide either --submission_basepath OR --input, not both."
    }

    if( !submission_basepath && !input_samplesheet ) {
        error "You must provide either --submission_basepath (a directory) or --input (a samplesheet CSV)."
    }

    if( submission_basepath ) {

    // first step: create samplesheet from metadata.json file
    METADATA_TO_SAMPLESHEET(
        Channel.fromPath(params.submission_basepath, checkIfExists: true)
    )

        ch_samplesheet = METADATA_TO_SAMPLESHEET.out.samplesheet
        ch_genome  = METADATA_TO_SAMPLESHEET.out.genome  
    
    } else if ( input_samplesheet ) {
        // Use the provided samplesheet file directly
        ch_samplesheet = Channel.of( params.input )
        ch_genome  = Channel.of ( params.genome )

    } 

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        ch_samplesheet
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_GRZQC (
        PIPELINE_INITIALISATION.out.samplesheet,
        ch_genome
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_GRZQC.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
