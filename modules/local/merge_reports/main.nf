process MERGE_REPORTS {

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed624a85396ad8cfe079da9b0bf12bf9822bbebcbbe926c24bb49906665ed4be/data'
        : 'community.wave.seqera.io/library/pip_gzip-utils_openpyxl_pandas:cd97ba68cc5b8463'}"

    input:
    path csv_files

    output:
    path "*report.csv"
    path "*report.xlsx"
    path "*report_mqc.csv", emit: multiqc
    path ('versions.yml'), emit: versions

    script:
    def prefix = task.ext.prefix ?: ""
    def create_alias = task.ext.create_alias ? true : false
    def alias = task.ext.create_alias ?: ''

    """
    merge_reports.py ${csv_files} --output_prefix ${prefix}report

    # If enabled, create a symbolic link named 'report.csv' pointing to the prefixed output
    if [ "${create_alias}" = true ]; then
        ln -s ${prefix}report.csv ${alias}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    def create_alias = task.ext.create_alias ? true : false
    def alias = task.ext.create_alias ?: ''
    """
    touch ${prefix}report.csv
    touch ${prefix}report.xlsx

    if [ "${create_alias}" = true ]; then
        ln -s ${prefix}report.csv ${alias}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
