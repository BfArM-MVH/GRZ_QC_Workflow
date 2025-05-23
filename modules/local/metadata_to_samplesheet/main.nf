process METADATA_TO_SAMPLESHEET {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/02/02bc9315ec85b45ba0576aaba63285c5908187067af3dcfd3a1ede24920cb8f9/data'
        : 'community.wave.seqera.io/library/grz-pydantic-models:1.4.0--570e8259614c00d0'}"

    input:
    path submission_basepath

    output:
    path ("*samplesheet.csv"), emit: samplesheet

    script:
    """
    metadata_to_samplesheet.py "${submission_basepath}"
    """
}
