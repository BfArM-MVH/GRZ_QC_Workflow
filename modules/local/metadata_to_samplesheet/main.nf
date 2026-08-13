process METADATA_TO_SAMPLESHEET {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/42/427bd6491745418d288c4bdfb366b63f83ff7211c5b39c79f5e21e1646acb8e7/data'
        : 'community.wave.seqera.io/library/pandas_grz-pydantic-models:2ab7f5e78743d861'}"

    input:
    path submission_basepath

    output:
    path ("*samplesheet.csv"), emit: samplesheet

    script:
    """
    metadata_to_samplesheet.py "${submission_basepath}"
    """
}
