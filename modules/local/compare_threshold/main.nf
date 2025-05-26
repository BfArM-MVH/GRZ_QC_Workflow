process COMPARE_THRESHOLD {

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed624a85396ad8cfe079da9b0bf12bf9822bbebcbbe926c24bb49906665ed4be/data'
        : 'community.wave.seqera.io/library/pip_gzip-utils_openpyxl_pandas:cd97ba68cc5b8463'}"

    input:
    tuple val(meta), path(summary), path(bed), path(fastp_jsons)
    path thresholds

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
        --fastp_json ${fastp_jsons} \\
        --mosdepth_global_summary ${summary} \\
        --mosdepth_target_regions_bed ${bed} \\
        --thresholds ${thresholds} \\
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
