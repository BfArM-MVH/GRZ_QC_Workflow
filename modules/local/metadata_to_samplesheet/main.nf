
process METADATA_TO_SAMPLESHEET {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ae/ae18ec651e9014c4d403a2837f348f4042b09ff565a50db2b01ba9b3344dc046/data' :
        'community.wave.seqera.io/library/gzip_jq_openpyxl_pandas:24ed4e1917ca1d2f' }"

    input:
    path submission_basepath

    output:
    path("*samplesheet.csv"), emit: samplesheet
    env("genome"),  emit: genome

    // when:
    // task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metadata_to_samplesheet.py "${submission_basepath}" 
    genome=\$(jq -r '.donors[0].labData[0].sequenceData.referenceGenome' "${submission_basepath}/metadata/metadata.json")
    """

}
