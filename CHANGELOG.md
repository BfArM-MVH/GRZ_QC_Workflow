# BfArM-MVH/GRZ_QC_Workflow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0-rc1 - [29.04.2025]

### `Added`
- Mark duplicates 
- decrease reference reading to one (either GRCh37 or GRCh38)
- two inputs: metadata path or samplesheet
- observed resource usage on README
- Save workflow reports from CI
- update README about usage of costum config + update base.config for bwamem

### `Fixed`
- multiqc reports of mosdepth and fastp
- change fractionBasesAboveQualityThreshold to percentBasesAboveQualityThreshold following meta json file
- picard TMP_DIR

### `Dependencies`

### `Deprecated`

## v0.2.0-alpha - [19.03.2025]

Initial release of nf-core/grzqc, created with the nf-core template.