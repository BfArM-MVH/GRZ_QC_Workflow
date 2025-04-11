[![GitHub Actions CI Status](https://github.com/nf-core/grzqc/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/grzqc/actions/workflows/ci.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**BfArM-MVH/GRZ_QC_Workflow** performs extended quality control of GRZ submissions according to the defined thresholds.

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [`FASTP`](https://github.com/OpenGene/fastp))
2. Alignment using ([`BWAMEM2`](https://github.com/bwa-mem2/bwa-mem2))
3. MarkDuplicates using ([`Picard`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard))
4. Coverage calculation by ([`Mosdepth`](https://github.com/brentp/mosdepth))
5. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

For the details on coverage calculation of different library types, see the [documentation](docs/details.md).

For the exact command lines executed by the pipeline you may check the workflow reports automatically generated by the test pipeline on GitHub.
To access these reports, click on the title of the latest successful [pipeline run](https://github.com/BfArM-MVH/GRZ_QC_Workflow/actions/workflows/main.yml) and download the `nextflow-pipeline-info` artifact under the "Artifacts" section at the bottom of the page.
The command lines are detailed the the "Tasks" table at the bottom of the `execution_report_*.html` file.

## Setup

- Install nextflow (and dependencies)
- Make sure to have either conda, docker or singularity.
- Clone the github repository

```bash
git clone https://github.com/BfArM-MVH/GRZ_QC_Workflow.git
output_basepath="path/to/analysis/dir"
mkdir -p ${output_basepath}/grzqc_output
```

This pipeline will automatically download the necessary _reference genomes_ and creates an _BWA index_ from them.
However, when running this pipeline multiple times on different submissions, the download and indexing steps create _unnecessary overhead_.

Therefore, it is recommended to run test_GRCh37 and test_GRCh38 profiles to set-up reference files and to be sure if all the necessary images and containers are settles-up correctly. Thus, depending on the reference genome you will use, run the profile:

```bash
nextflow run main.nf \
    -profile test_GRCh37,docker
    --save_reference_path "${output_basepath}"
```

or/and

```bash
nextflow run main.nf \
    -profile test_GRCh38,docker
    --save_reference_path "${output_basepath}"
```

\*\*Please use replace docker with singularity or conda depending on your system. This pipeline is able to run all profiles.

Now, all the necessary files are saved into _output_/_basepath_/references. In the next section you will see how can use these files and avoid rerunning the genome downloading and indexing steps.

and you can delete test results safely:

```bash
rm -rf ${projectDir}/tests/results
```

## Usage
This pipeline needs one of the following two inputs:

1. a submission base directory path with a folder structure following [GRZ submission standard](https://github.com/BfArM-MVH/grz-cli?tab=readme-ov-file#introduction). You can also check [test datasets](https://www.cmm.in.tum.de/public/grz-example-submissions/).

You can run the pipeline using:

```bash
submission_basepath="path/to/submission/base/directory"
nextflow run main.nf \
    -profile docker \
    --outdir "${output_basepath}/grzqc_output/" \
    --submission_basepath "${submission_basepath}" \
    -c conf/grzqc_GRCh37.config
```

or with `-c conf/grzqc_GRCh38.config` flag for GRCh38.

If you copy the _reference genomes_ and _BWA index_ somewhere else after the test run, you can also change the lines in `conf/grzqc_GRCh37.config` and `conf/grzqc_GRCh38.config`.

```bash
    fasta = "${outdir}/../reference/GRCh37/genome.fasta"
    fai   = "${outdir}/../reference/GRCh37/genome.fasta.fai"
    bwa   = "${outdir}/../reference/GRCh37/bwamem2"
```

2. In addition to use a submission base directory path, you can also use a csv samplesheet as input. This gives more flexibilty, as you don't need a GRZ submission directory. See the [documentation](docs/usage.md).

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

For more details about the output files and reports, please refer to the [output documentation](docs/output.md).

## Thresholds

QC thresholds are read from [`thresholds.json`](https://github.com/BfArM-MVH/GRZ_QC_Workflow/blob/main/assets/default_files/thresholds.json), which uses the values [defined by BfArM](https://www.bfarm.de/SharedDocs/Downloads/DE/Forschung/modellvorhaben-genomsequenzierung/Qs-durch-GRZ.pdf?__blob=publicationFile).

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

The default reference files used in this pipeline is part of AWS iGenomes. Please follow the instructions [here](https://ewels.github.io/AWS-iGenomes/) to download:

- s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/
- s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/

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

**BfArM-MVH/GRZ_QC_Workflow** was originally written by Yun Wang, Kübra Narci, Shounak Chakraborty and Florian R. Hölzlwimmer.
