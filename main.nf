#!/usr/bin/env nextflow

include { GRZQC                   } from './workflows/grzqc'
include { paramsSummaryLog        } from 'plugin/nf-schema'
include { validateParameters      } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { METADATA_TO_SAMPLESHEET } from './modules/local/metadata_to_samplesheet'

workflow {
    //
    // Check whether a submission_basepath or a samplesheet is provided
    //
    def submission_basepath = params.submission_basepath && params.submission_basepath.trim() != 'null' ? params.submission_basepath : null
    def input_samplesheet = params.input && params.input.trim() != 'null' ? params.input : null

    if (submission_basepath && input_samplesheet) {
        error("Please provide either --submission_basepath OR --input, not both.")
    }

    if (!submission_basepath && !input_samplesheet) {
        error("You must provide either --submission_basepath (a directory) or --input (a samplesheet CSV).")
    }

    if (submission_basepath) {

        // first step: create samplesheet from metadata.json file
        METADATA_TO_SAMPLESHEET(
            Channel.fromPath(params.submission_basepath, checkIfExists: true)
        )

        ch_samplesheet = METADATA_TO_SAMPLESHEET.out.samplesheet
    }
    else if (input_samplesheet) {
        // Use the provided samplesheet file directly
        ch_samplesheet = Channel.of(params.input)
    }

    //
    // Print parameter summary to stdout. This will display the parameters
    // that differ from the default given in the JSON schema
    //
    log.info(paramsSummaryLog(workflow))

    //
    // Validate the parameters using nextflow_schema.json or the schema
    // given via the validation.parametersSchema configuration option
    //
    if (params.validate_params) {
        validateParameters()
    }

    // construct the samplesheet channel using the schema
    // read_group is necessary for BWAMEM2_MEM
    ch_samplesheet
        .flatMap { samplesheet ->
            samplesheetToList(samplesheet.toString(), "${projectDir}/assets/schema_input.json")
        }
        .map { meta, reads1, reads2, alignment ->
            def id = meta.runId ? "${meta.sample}_${meta.runId}" : "${meta.sample}"
            def read_group = "@RG\\tID:${id}\\tPL:${meta.sequencer}\\tSM:${meta.sample}"
            def single_end = !reads2
            def reads = single_end ? [reads1] : [reads1, reads2]

            return [
                meta + [
                    id: id,
                    single_end: single_end,
                    read_group: read_group,
                ],
                alignment ? [] : reads,
                alignment,
            ]
        }
        .set { ch_samplesheet }

    if (params.genome) {
        ch_genome = Channel.of(params.genome)
    }
    else {
        ch_genome = ch_samplesheet.map { meta, _reads, _alignments -> meta.reference }.unique().first()
    }

    GRZQC(
        ch_samplesheet,
        ch_genome,
    )
}
