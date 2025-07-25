process METADATA_TO_SAMPLESHEET {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/93adbb05c1b65c9c4fdc43b5142ea9a045a0373c37e1ec1097d47c8d8c01c405/data'
        : 'community.wave.seqera.io/library/grz-pydantic-models_pandas:d192123ee8af7fc4'}"

    input:
    path submission_basepath

    output:
    path ("*samplesheet.csv"), emit: samplesheet

    script:
    """
    metadata_to_samplesheet.py "${submission_basepath}"
    """
}
