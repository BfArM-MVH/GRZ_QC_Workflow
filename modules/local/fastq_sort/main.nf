process FASTQ_SORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6b/6b85d6776227c966f456dc4b1d1edf6db645bef4cfce2bdc928151f64fca5316/data'
        : 'community.wave.seqera.io/library/fastq-tools_samtools_parallel:2a78b0b4ab3491b2'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.sorted.fastq.gz'), emit: reads
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cpus_per_subtask = Math.max(task.cpus.intdiv(reads.size()), 1)

    if (meta.single_end) {
        // single-end FASTQs don't require sorting
        """
        ln -s ${reads} ${prefix}.sorted.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            parallel: \$(parallel --version | head -n1 | sed 's/^GNU parallel //')
            fastq-sort: \$(fastq-sort --version | cut -f3 -d' ')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
    else {
        """
        parallel 'bgzip --stdout --decompress --threads ${cpus_per_subtask} {} | fastq-sort ${args} | bgzip --stdout --threads ${cpus_per_subtask} > ${prefix}.R{#}.sorted.fastq.gz' ::: ${reads.join(' ')}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            parallel: \$(parallel --version | head -n1 | sed 's/^GNU parallel //')
            fastq-sort: \$(fastq-sort --version | cut -f3 -d' ')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
