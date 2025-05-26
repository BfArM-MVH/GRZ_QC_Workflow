process METADATA_TO_SAMPLESHEET {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a9/a9edc101b57b3042094a569950f2d0772333156ef83984bb24d701d71bd22030/data'
        : 'community.wave.seqera.io/library/grz-pydantic-models_pandas:b31707c70fa2229e'}"

    input:
    path submission_basepath

    output:
    path ("*samplesheet.csv"), emit: samplesheet

    script:
    """
    metadata_to_samplesheet.py "${submission_basepath}"
    """
}
