process CONVERT_BED_CHROM {

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed624a85396ad8cfe079da9b0bf12bf9822bbebcbbe926c24bb49906665ed4be/data'
        : 'community.wave.seqera.io/library/pip_gzip-utils_openpyxl_pandas:cd97ba68cc5b8463'}"

    input:
    tuple val(meta), path(bed_file)
    path mapping_file

    output:
    tuple val(meta), path("*.converted_bed.bed"), emit: converted_bed
    path ('versions.yml'), emit: versions

    script:
    """
    convert_bed_chrom.py \\
        --bed ${bed_file} \\
        --mapping ${mapping_file} \\
        --output "${meta.id}.converted_bed.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.converted_bed.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
