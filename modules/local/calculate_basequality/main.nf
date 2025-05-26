process CALCULATE_BASEQUALITY {

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ba/ba5dc444f779fd6c19aaf6b5e767e441f4578b3d8841280a5b6e6a3a23a05044/data'
        : 'ccommunity.wave.seqera.io/library/pip_pysam:0756c7708ed1598b'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}.json"), emit: json
    path ('versions.yml'), emit: versions

    script:
    """
    #!/usr/bin/env bash
    if [ -f "${meta.fastp_json}" ]; then
        cp "${meta.fastp_json}" "${meta.id}.json"
    else
        calculate_basequality.py \\
            --input "${input}" \\
            --output "${meta.id}.json"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
