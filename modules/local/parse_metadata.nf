process PARSE_METADATA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975/data' :
        'community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264' }"

    input:
    tuple val(meta), path(submission_base_path, stageAs: "input*/*")

    output:
    tuple val(meta), path("*samplesheet.csv"), emit: samplesheet

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metadata_to_samplesheet.py $submission_base_path
    """

}
