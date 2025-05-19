#!/usr/bin/env nextflow

include { GRZQC                   } from './workflows/grzqc'
include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { METADATA_TO_SAMPLESHEET } from './modules/local/metadata_to_samplesheet'

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
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //
    log.info paramsSummaryLog(workflow)

    //
    // Validate the parameters using nextflow_schema.json or the schema
    // given via the validation.parametersSchema configuration option
    //
    if(params.validate_params) {
      validateParameters()
    }

    // construct the samplesheet channel using the schema
    ch_samplesheet
        .flatMap  { 
            samplesheet ->
                samplesheetToList(samplesheet.toString(), "${projectDir}/assets/schema_input.json")
        }
        .map {
            meta, fastq_1, fastq_2, bed_file ->
                if (!fastq_2) {
                    return [ meta + [ single_end:true ], [ fastq_1 ], bed_file]
                } else {
                    return [ meta + [ single_end:false ], [ fastq_1, fastq_2 ] , bed_file]
                }
        }
        .groupTuple()
        .map {
            meta, fastqs, bed_file ->
                return [ meta, fastqs.flatten(), bed_file ]
        }
        .set{ ch_samplesheet }

    GRZQC (
        ch_samplesheet,
        ch_genome
    )
}
