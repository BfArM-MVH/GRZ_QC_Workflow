[![GitHub Actions CI Status](https://github.com/nf-core/grzqc/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/grzqc/actions/workflows/ci.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**BfArM-MVH/GRZ_QC_Workflow** performs extended quality control of GRZ submissions according to the defined thresholds.

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [`FASTP`](https://github.com/OpenGene/fastp))
2. Alignment using ([`BWAMEM2`](https://github.com/bwa-mem2/bwa-mem2))
3. Coverage calculation by ([`Mosdepth`](https://github.com/brentp/mosdepth))
4. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Setup

- Install nextflow (and dependencies)
- Make sure to have either conda, docker or singularity.
- Clone the github repository

```bash
git clone https://github.com/BfArM-MVH/GRZ_QC_Workflow.git
$output_path = "path/to/analysis/dir"
```

### Setting up reference files

This pipeline will automatically download the necessary reference genomes and creates an BWA index from them.
However, when running this pipeline multiple times on different submissions, the download and indexing steps create unnecessary overhead.

To skip downloading the reference genomes, you can also download the necessary reference genome FASTA files to some shared location:

You can download files here

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz $shared_directory/references
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
mv hg38.fa.gz $shared_directory/references
```

or from (igenomes)[https://ewels.github.io/AWS-iGenomes/]: Homo sapiens : UCSC : hg38/hg19 : FASTA

Then you can update the file paths in `conf/grzqc.conf`:

```bash
params {
    [...]
    fasta = "$shared_directory/references/hg19.fa.gz"
    or
    fasta = "$shared_directory/references/hg38.fa.gz"
}
```

by replacing `$shared_directory` with the absolute path to the shared directory.

After the first run, you can also copy the BWAMEM2 index to the shared directory:

```bash
cp -r "${output_basepath}/grzqc_output/references/" "$shared_directory/references/"
```

or from (igenomes)[https://ewels.github.io/AWS-iGenomes/]: Homo sapiens : UCSC : hg38/hg19 : BWA

and configure it in `conf/grzqc.conf`:

```bash
params {
    [...]
    bwa = "$shared_directory/references/GRCh37/bwamem2"
    or
    bwa = "$shared_directory/references/GRCh38/bwamem2"

}
```

by replacing `$shared_directory` with the absolute path to the shared directory.

## Usage

This pipeline needs a samplesheet which is generated automatically from the metadata.json file included in the submission base directory. Please make sure that the submission base directory has the required folder structure. The script `run_grzqc.sh` parses the metadata.json file to create a nextflow samplesheet:

```bash
python3 bin/metadata_to_samplesheet.py \
    "${submission_basepath}" \
    "${output_basepath}/grzqc_output/grzqc_samplesheet.csv"
```

Now, you can run the pipeline using:

```bash
nextflow run main.nf \
    -profile docker \
    --outdir "${output_basepath}/grzqc_output/" \
    --input "${output_basepath}/grzqc_output/grzqc_samplesheet.csv" \
    --genome "hg38"

```

For your next run, you can use prebuild references. Please prepare your own config file to do so.

## Pipeline output

Output :

| Column                                       | Description                                                             |
| -------------------------------------------- | ----------------------------------------------------------------------- |
| `sampleId`                                   | Sample ID                                                               |
| `labDataName`                                | Lab data name                                                           |
| `libraryType`                                | Library type, e.g., `wes` for whole-exome sequencing                    |
| `sequenceSubtype`                            | Sequence subtype, e.g., `somatic` or `germline`                         |
| `genomicStudySubtype`                        | Genomic study subtype, e.g., `tumor+germline`                           |
| `meanDepthOfCoverage`                        | Mean depth of coverage                                                  |
| `meanDepthOfCoverageRequired`                | Mean depth of coverage required to pass QC                              |
| `fractionBasesAboveQualityThreshold`         | Fraction of bases passing the quality threshold                         |
| `qualityThreshold`                           | The quality threshold to pass                                           |
| `fractionBasesAboveQualityThresholdRequired` | Fraction of bases above the quality threshold required to pass QC       |
| `targetedRegionsAboveMinCoverage`            | Fraction of targeted regions above minimum coverage                     |
| `minCoverage`                                | Minimum coverage for target regions                                     |
| `targetedRegionsAboveMinCoverageRequired`    | Fraction of targeted regions above minimum coverage required to pass QC |
| `passedQC`                                   | `true` when QC passed, otherwise `false`                                |

## Running the pipeline offline

Nextflow can automatically retrieve almost everything necessary to execute a pipeline from the web, including pipeline code, software dependencies, reference genomes, and remote data sources.

However, if your analysis must run on a system without _internet access_, you will need to take a few additional steps to ensure all required components are available locally. First, download everything on an internet-connected system (such as your personal computer) and then transfer the files to the offline system using your preferred method.

To set up an offline environment, you will need three key components: a functioning Nextflow installation, the pipeline assets, and the required reference genomes.

On a computer with an internet connection, to download the pipeline, run:

```bash
nf-core pipelines download BfArM-MVH/GRZ_QC_Workflow
```

Add the argument `--container-system singularity` to also fetch the singularity container(s).

Then download the necessary plugins and lace it under `${NXF_HOME}/plugins`:

```bash
nextflow plugin install nf-schema@2.1.1

```

For more detailed information please check ["Running offline by nf-core"](https://nf-co.re/docs/usage/getting_started/offline)


## Estimated resource requirements

Using the 466 GB `WGS_tumor+germline` test submission dataset from the
[example GRZ submissions](https://www.cmm.in.tum.de/public/grz-example-submissions),
the pipeline used the following resources:

- 828 CPU hours
- 72 GB maximum RAM (genome indexing)
- 2 TB storage (including the input files)

The biggest jobs were the two bwa-mem2 alignments which used 300 CPU hours each
and a maximum of 48 GB of RAM.

## Contributions and Support

**BfArM-MVH/GRZ_QC_Workflow** was originally written by Shounak Chakraborty, Yun Wang, Kübra Narci and Florian R. Hölzlwimmer.
