process SAVE_REFERENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'quay.io/nf-core/ubuntu:20.04'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(index)
    val genome

    output:
    tuple val(meta), path("reference/${genome}"), emit: reference

    script:
    """
    #!/bin/bash

    mkdir -p reference/${genome}
    cp ${fasta} reference/${genome}
    cp ${fai} reference/${genome}
    cp -r ${index} reference/${genome}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    """
    #!/bin/bash
    mkdir reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
