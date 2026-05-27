process METADATA_TO_SAMPLESHEET {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/97/97e40ef01757e3804562c0e7de9ada44dd8c4b6cdcee3e4d5d554aca68b22545/data'
        : 'community.wave.seqera.io/library/pandas_grz-pydantic-models:1e97fd2149afc157'}"

    input:
    path submission_basepath

    output:
    path ("*samplesheet.csv"), emit: samplesheet

    script:
    """
    metadata_to_samplesheet.py "${submission_basepath}"
    """
}
