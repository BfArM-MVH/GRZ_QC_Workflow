
process MERGE_REPORTS {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed624a85396ad8cfe079da9b0bf12bf9822bbebcbbe926c24bb49906665ed4be/data' :
        'community.wave.seqera.io/library/pip_gzip-utils_openpyxl_pandas:b1a85fb63f75244d' }"

    input:
        path (csv_files)

    output:
        path "merged_result.csv"
        path "merged_result.xlsx"
        path('versions.yml')      , emit: versions

    script:
    """
    merge_reports.py $csv_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    stub:
    """
    touch merged_result.csv
    touch merged_result.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
