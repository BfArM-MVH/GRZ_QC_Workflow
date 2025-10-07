process METADATA_TO_SAMPLESHEET {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1a/1abe0e2c1653fe6b64fcfd9b3ac03a3572897609a745f67f6b306e73207f1c81/data'
        : 'community.wave.seqera.io/library/grz-pydantic-models_pandas:0027fb9247a1e699'}"

    input:
    path submission_basepath

    output:
    path ("*samplesheet.csv"), emit: samplesheet

    script:
    """
    metadata_to_samplesheet.py "${submission_basepath}"
    """
}
