process MERGE_REPORTS {

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed624a85396ad8cfe079da9b0bf12bf9822bbebcbbe926c24bb49906665ed4be/data'
        : 'community.wave.seqera.io/library/pip_gzip-utils_openpyxl_pandas:cd97ba68cc5b8463'}"

    input:
    path csv_files

    output:
    path "report.csv"
    path "report.xlsx"
    path ('versions.yml'), emit: versions

    script:
    """
    merge_reports.py ${csv_files} --output_prefix "report"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch report.csv
    touch report.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
