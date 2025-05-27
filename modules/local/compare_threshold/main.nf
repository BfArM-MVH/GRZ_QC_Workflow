process COMPARE_THRESHOLD {

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a9/a9edc101b57b3042094a569950f2d0772333156ef83984bb24d701d71bd22030/data'
        : 'community.wave.seqera.io/library/grz-pydantic-models_pandas:b31707c70fa2229e'}"

    input:
    tuple val(meta), path(summary), path(bed), path(fastp_jsons)

    output:
    path ('*.result.csv'), emit: result_csv
    path ('versions.yml'), emit: versions

    script:
    """
    compare_threshold.py \\
        --sample_id ${meta.id} \\
        --labDataName "${meta.labDataName}" \\
        --libraryType "${meta.libraryType}" \\
        --sequenceSubtype "${meta.sequenceSubtype}" \\
        --genomicStudySubtype "${meta.genomicStudySubtype}" \\
        --meanDepthOfCoverage "${meta.meanDepthOfCoverage}" \\
        --targetedRegionsAboveMinCoverage "${meta.targetedRegionsAboveMinCoverage}" \\
        --percentBasesAboveQualityThreshold "${meta.percentBasesAboveQualityThreshold}" \\
        --fastp_json ${fastp_jsons} \\
        --mosdepth_global_summary ${summary} \\
        --mosdepth_target_regions_bed ${bed} \\
        --output ${meta.id}.result.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.result.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
