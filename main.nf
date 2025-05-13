include { PARSE_SUBMISSION } from './modules/local/submission/parse'
include { CAT_FASTQ as CAT_FASTQ_SRT; CAT_FASTQ as CAT_FASTQ_LNG } from './modules/nf-core/cat/fastq'
include { FASTQC } from './modules/nf-core/fastqc/main'
include { FASTP } from './modules/nf-core/fastp/main'
include { BWAMEM2_INDEX } from './modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM } from './modules/nf-core/bwamem2/mem/main'
include { PICARD_ADDORREPLACEREADGROUPS } from './modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard/markduplicates/main'
include { PBTK_PBINDEX } from './modules/nf-core/pbtk/pbindex/main'
include { PBTK_BAM2FASTQ } from './modules/nf-core/pbtk/bam2fastq/main'
include { FASTPLONG } from './modules/local/fastplong'
include { MINIMAP2_INDEX } from './modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN } from './modules/nf-core/minimap2/align/main'
include { MOSDEPTH } from './modules/nf-core/mosdepth/main'
include { CHECK_THRESHOLDS } from './modules/local/thresholds'
include { MULTIQC } from './modules/nf-core/multiqc/main'


workflow {
  def submission_root = file(params.submission)
  def samples = PARSE_SUBMISSION(submission_root)
    .samples
    .map{ sample ->
        sample.regions = sample.regions == null ? file(params.genomes[sample.meta.reference_genome].regions) : sample.regions
        sample
    }
    .dump(tag: 'samples', pretty: true)
    .branch {
      lng: it.meta.library_type.endsWith('_lr')
      srt: true
    }
  def genomes = channel.of(params.genomes)
    .flatMap{ it.entrySet() }
    .map{ [it.getKey(), it.getValue()] }

  // short read workflow
  BWAMEM2_INDEX(
    samples.srt.map{ [it.meta.reference_genome, it.meta] }
      .join( genomes )
      .map{ [["id": it[0]], file(it[2].fasta)] }
  )

  samples.srt
    .map{ [it.meta, it.reads] }
    | CAT_FASTQ_SRT

  FASTQC(CAT_FASTQ_SRT.out.reads)

  FASTP(
    CAT_FASTQ_SRT.out.reads,
    [],    // adapter sequences
    false, // don't discard passed reads
    false, // don't save failed reads
    false  // don't save reads to merged file
  )

  BWAMEM2_MEM(
    FASTP.out.reads,
    BWAMEM2_INDEX.out.index
      .map{ [it[0].id, it[1]] }
      .combine( FASTP.out.reads.map{ [it[0].reference_genome, it] }, by: 0)
      .map{ [it[2][0], it[1]] },
    genomes.map{ [it[0], file(it[1].fasta)] }
      .combine( FASTP.out.reads.map{ [it[0].reference_genome, it] }, by: 0)
      .map{ [it[2][0], it[1]] },
    true  // sort alignments
  )

  PICARD_ADDORREPLACEREADGROUPS(
    BWAMEM2_MEM.out.bam,
    channel.value([[:], []]),  // no reference needed for BAM input
    channel.value([[:], []])   // ^
  )

  PICARD_MARKDUPLICATES(
    PICARD_ADDORREPLACEREADGROUPS.out.bam,
    channel.value([[:], []]),  // no reference needed for BAM input
    channel.value([[:], []])   // ^
  )

  // long read workflow
  def subsamples_lng = samples.lng.multiMap { sample ->
    fastq: [
      sample.meta,
      sample.reads
        .withIndex()
        .collect{ [sample.read_file_types[it[1]] == "fastq", it[0]] }
        .grep{ it[0] }
        .collect{ it[1] }
    ]
    bam: [
      sample.meta,
      sample.reads
        .withIndex()
        .collect{ [sample.read_file_types[it[1]] == "bam", it[0]] }
        .grep{ it[0] }
        .collect{ it[1] }
    ]
  }

  def bams = subsamples_lng.bam.filter{ it[1].size() }.flatMap { it[1].collect{ v -> [["id": it[0].id + "_" + v.name, "sample": it[0].id], v] } }

  PBTK_PBINDEX(bams)

  PBTK_BAM2FASTQ(
    bams.map{ [it[0].id, it] }.join(
      PBTK_PBINDEX.out.pbi.map{ [it[0].id, it] }
    ).map{ [it[1][0], it[1][1], it[2][1]] }
  )

  CAT_FASTQ_LNG(
    subsamples_lng.fastq.map{ [it[0].id, it] }.join(
      PBTK_BAM2FASTQ.out.fastq.map{ [it[0].sample, it[1]] }.groupTuple(),
      remainder: true
    ).map{ [it[1][0], (it[1][1] + it[2]).grep{ it != null }] }
  )

  FASTPLONG(
    CAT_FASTQ_LNG.out.reads,
    channel.value([]),  // optional adapters
    false,              // don't discard passing reads
    false               // don't save failed reads
  )

  MINIMAP2_INDEX(
    FASTPLONG.out.reads.map{ [it[0].reference_genome] }
      .join( genomes )
      .map{ [["id": it[0]], file(it[1].fasta)] }
  )

  MINIMAP2_ALIGN(
    FASTPLONG.out.reads,
    MINIMAP2_INDEX.out.index
      .map{ [it[0].id, it[1]] }
      .combine( FASTPLONG.out.reads.map{ [it[0].reference_genome, it] }, by: 0 )
      .map{ [it[2][0], it[1]] },
    true,   // output BAM
    "bai",  // BAM index extension
    false,  // PAF CIGAR output, ignored with BAM output
    true    // output CIGAR in CG tag in BAM
  )

  // combined read workflow
  def bams_srt = PICARD_MARKDUPLICATES.out.bam.map{ [it[0].id, it[0], it[1]] }.join(
    PICARD_MARKDUPLICATES.out.bai.map{ [it[0].id, it[1]] }
  )
  def bams_lng = MINIMAP2_ALIGN.out.bam.map{ [it[0].id, it[0], it[1]] }.join(
    MINIMAP2_ALIGN.out.index.map{ [it[0].id, it[1]] }
  )

  MOSDEPTH(
    bams_srt.mix(bams_lng)
      .join(samples.srt.mix(samples.lng).map{ sample -> [sample.meta.id, sample.regions] })
      .map{ it.subList(1, 5) },
    channel.value([[:], []])   // no reference needed for BAM input
  )

  CHECK_THRESHOLDS(
    samples.srt.mix(samples.lng).map{ sample -> [sample.meta.id, sample.meta] }
      .join(FASTP.out.json.mix(FASTPLONG.out.json).map{ meta, json -> [meta.id, json] })
      .join(MOSDEPTH.out.summary_txt.map{ meta, txt -> [meta.id, txt] })
      .join(MOSDEPTH.out.regions_bed.map{ meta, bed -> [meta.id, bed] })
      .map{ it.subList(1, 5) }
  )

  MULTIQC(
    FASTQC.out.zip.collect{ meta, zip -> zip }
      .mix(FASTP.out.json.collect{ meta, json -> json })
      .mix(FASTPLONG.out.json.collect{ meta, json -> json })
      .mix(MOSDEPTH.out.summary_txt.collect{ meta, txt -> txt })
      .mix(MOSDEPTH.out.global_txt.collect{ meta, txt -> txt })
      .mix(MOSDEPTH.out.regions_txt.collect{ meta, txt -> txt })
      .mix(CHECK_THRESHOLDS.out)
      .collect(),
    [],  // optional config
    [],  // optional extra config
    [],  // optional logo
    [],  // optional sample name replacement file
    []   // optional sample names
  )
}
