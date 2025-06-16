include { SAMTOOLS_BGZIP } from '../../../modules/nf-core/samtools/bgzip/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx'
include { BWAMEM2_INDEX  } from '../../../modules/nf-core/bwamem2/index'
include { MINIMAP2_INDEX } from '../../../modules/nf-core/minimap2/index'

workflow PREPARE_REFERENCES {
    take:
    ch_samplesheet_reads_srt
    ch_samplesheet_reads_lng
    ch_genome

    main:
    // general order of precedence: reference_path -> param -> default
    def reference_path = params.reference_path
    def ch_versions = Channel.empty()

    // prepare fasta
    def defaultFasta = [
        "GRCh37": "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa",
        "GRCh38": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz",
    ]

    def ch_fasta_state = ch_genome.branch { id ->
        cached: reference_path ? fasta_cached(id, reference_path).exists() : false
        not_cached: true
    }
    def ch_fasta_cached = ch_fasta_state.cached.map { id -> [[id: id], fasta_cached(id, reference_path)] }
    def ch_fasta_not_cached = ch_fasta_state.not_cached.map { id ->
        [[id: id], file(params.fasta ?: defaultFasta[id])]
    }
    SAMTOOLS_BGZIP(ch_fasta_not_cached)
    ch_versions = ch_versions.mix(SAMTOOLS_BGZIP.out.versions)
    def ch_fasta = ch_fasta_cached.mix(SAMTOOLS_BGZIP.out.fasta).dump(tag: 'fasta')

    // prepare fai
    def ch_fai_state = ch_fasta.branch { meta, _fasta ->
        cached: reference_path ? fai_cached(meta.id, reference_path).exists() : false
        not_cached: true
    }
    def ch_fai_cached = ch_fai_state.cached.map { meta, _fasta -> [meta, fai_cached(meta.id, reference_path)] }
    def faidx_fai = [[], []]
    def faidx_get_sizes = false
    SAMTOOLS_FAIDX(
        ch_fai_state.not_cached.map { meta, _fai -> [meta] }.join(ch_fasta),
        faidx_fai,
        faidx_get_sizes,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    def ch_fai = ch_fai_cached.mix(SAMTOOLS_FAIDX.out.fai).dump(tag: 'fai')

    // prepare bwa
    def ch_bwa_state = ch_fasta.branch { meta, _fasta ->
        cached: reference_path ? bwa_cached(meta.id, reference_path).exists() : false
        not_cached: true
    }
    def ch_bwa_cached = ch_bwa_state.cached.map { meta, _fasta -> [meta, bwa_cached(meta.id, reference_path)] }
    // only index with bwamem2 if there are short reads in samplesheet
    def ch_fasta_if_have_short_reads = ch_fasta.combine(ch_samplesheet_reads_srt).map { meta_fasta, path_fasta, _meta_row, _reads -> [meta_fasta, path_fasta] }.first()
    BWAMEM2_INDEX(
        ch_bwa_state.not_cached.map { meta, _bwa -> [meta] }.join(ch_fasta_if_have_short_reads)
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    def ch_bwa = ch_bwa_cached.mix(BWAMEM2_INDEX.out.index).dump(tag: 'bwa')

    // prepare mmi
    def ch_mmi_state = ch_fasta.branch { meta, _fasta ->
        cached: reference_path ? mmi_cached(meta.id, reference_path).exists() : false
        not_cached: true
    }
    def ch_mmi_cached = ch_mmi_state.cached.map { meta, _fasta -> [meta, mmi_cached(meta.id, reference_path)] }
    // only index with minimap2 if there are long reads in samplesheet
    def ch_fasta_if_have_long_reads = ch_fasta.combine(ch_samplesheet_reads_lng).map { meta_fasta, path_fasta, _meta_row, _reads -> [meta_fasta, path_fasta] }.first()
    MINIMAP2_INDEX(
        ch_mmi_state.not_cached.map { meta, _mmi -> [meta] }.join(ch_fasta_if_have_long_reads)
    )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    def ch_mmi = ch_mmi_cached.mix(MINIMAP2_INDEX.out.index).dump(tag: 'mmi')

    emit:
    fasta    = ch_fasta
    fai      = ch_fai
    bwa      = ch_bwa
    mmi      = ch_mmi
    versions = ch_versions
}

def fasta_cached(id, reference_path) {
    // get path that FASTA is cached to
    file(reference_path) / id / "${id}.fasta.gz"
}

def fai_cached(id, reference_path) {
    // get path that FAI is cached to
    file(reference_path) / id / "${id}.fasta.gz.fai"
}

def bwa_cached(id, reference_path) {
    // get path bwamem2 index is cached to
    file(reference_path) / id / "bwamem2/"
}

def mmi_cached(id, reference_path) {
    // get path minimap2 index is cached to
    file(reference_path) / id / "${id}.fasta.mmi"
}
