process COMPARE_THRESHOLD {
    tag "${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/135379b52c2d54842b471b5820082807fa0aae33cf1df118ebb3813dfe062c97/data'
        : 'community.wave.seqera.io/library/grz-pydantic-models_pandas:9c55bee92ebacc5d'}"

    input:
    tuple val(meta), path(summary), path(bed), path(fastp_jsons)

    output:
    tuple val(meta), path('*.result.csv'), emit: result_csv
    path ('versions.yml'), emit: versions

    script:
    def arg_meanDepthOfCoverageRequired = meta.meanDepthOfCoverageRequired ? "${meta.meanDepthOfCoverageRequired}" : '0'
    def arg_qualityThreshold = "${meta.qualityThreshold}"
    def arg_percentBasesAboveQualityThresholdRequired = meta.percentBasesAboveQualityThresholdRequired ? "${meta.percentBasesAboveQualityThresholdRequired}" : '0'
    def arg_minCoverage = meta.minCoverage ? "${meta.minCoverage}" : '0'
    def arg_targetedRegionsAboveMinCoverageRequired = meta.targetedRegionsAboveMinCoverageRequired ? "${meta.targetedRegionsAboveMinCoverageRequired}" : '0'

    def arg_meanDepthOfCoverage = meta.meanDepthOfCoverage ? "--meanDepthOfCoverage ${meta.meanDepthOfCoverage}" : ''
    def arg_targetedRegionsAboveMinCoverage = meta.targetedRegionsAboveMinCoverage ? "--targetedRegionsAboveMinCoverage ${meta.targetedRegionsAboveMinCoverage}" : ''
    def arg_percentBasesAboveQualityThreshold = meta.percentBasesAboveQualityThreshold ? "--percentBasesAboveQualityThreshold ${meta.percentBasesAboveQualityThreshold}" : ''
    def dedup = meta.is_deduplicated ? 'nonredundant' : 'redundant'
    // if using deduplicated data name output to redundant

    """
    compare_threshold.py \\
        --sample_id ${meta.id} \\
        --labDataName "${meta.labDataName ?: ''}" \\
        --donorPseudonym "${meta.donorPseudonym ?: ''}" \\
        --libraryType "${meta.libraryType}" \\
        --sequenceSubtype "${meta.sequenceSubtype}" \\
        --genomicStudySubtype "${meta.genomicStudySubtype}" \\
        --meanDepthOfCoverageRequired ${arg_meanDepthOfCoverageRequired} \\
        --qualityThreshold ${arg_qualityThreshold} \\
        --percentBasesAboveQualityThresholdRequired ${arg_percentBasesAboveQualityThresholdRequired} \\
        --minCoverage ${arg_minCoverage} \\
        --targetedRegionsAboveMinCoverageRequired ${arg_targetedRegionsAboveMinCoverageRequired} \\
        ${arg_meanDepthOfCoverage} \\
        ${arg_targetedRegionsAboveMinCoverage} \\
        ${arg_percentBasesAboveQualityThreshold} \\
        --fastp_json ${fastp_jsons} \\
        --mosdepth_global_summary ${summary} \\
        --mosdepth_target_regions_bed ${bed} \\
        --output ${meta.id}.${dedup}.result.csv

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
