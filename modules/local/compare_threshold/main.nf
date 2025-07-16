process COMPARE_THRESHOLD {
    tag "${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b249411a56a717bd02e6029ed0c540f170fe9d197dfa6ded0b2ce3581270c39/data'
        : 'community.wave.seqera.io/library/grz-pydantic-models_pandas:b9466bd2ba15f092'}"

    input:
    tuple val(meta), path(summary), path(bed), path(fastp_jsons)

    output:
    path ('*.result.csv'), emit: result_csv
    path ('versions.yml'), emit: versions

    script:
    def arg_meanDepthOfCoverage = meta.meanDepthOfCoverage ? "--meanDepthOfCoverage ${meta.meanDepthOfCoverage}" : ''
    def arg_targetedRegionsAboveMinCoverage = meta.targetedRegionsAboveMinCoverage ? "--targetedRegionsAboveMinCoverage ${meta.targetedRegionsAboveMinCoverage}" : ''
    def arg_percentBasesAboveQualityThreshold = meta.percentBasesAboveQualityThreshold ? "--percentBasesAboveQualityThreshold ${meta.percentBasesAboveQualityThreshold}" : ''

    """
    compare_threshold.py \\
        --sample_id ${meta.id} \\
        --labDataName "${meta.labDataName ?: ''}" \\
        --donorPseudonym "${meta.donorPseudonym ?: ''}" \\
        --libraryType "${meta.libraryType}" \\
        --sequenceSubtype "${meta.sequenceSubtype}" \\
        --genomicStudySubtype "${meta.genomicStudySubtype}" \\
        ${arg_meanDepthOfCoverage} \\
        ${arg_targetedRegionsAboveMinCoverage} \\
        ${arg_percentBasesAboveQualityThreshold} \\
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
