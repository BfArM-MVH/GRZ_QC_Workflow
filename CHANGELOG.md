# BfArM-MVH/GRZ_QC_Workflow: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0-rc1 - [29.04.2025]

### `Added`

- Mark duplicates [#41](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/41)
- decrease reference reading to one (either GRCh37 or GRCh38) [#36](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/36)
- two inputs: metadata path or samplesheet [#54](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/54)
- observed resource usage on README [#35](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/35)
- Save workflow reports from CI [#48](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/48)
- update README about usage of costum config + update base.config for bwamem [#63](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/63)
- Add disclaimer clarifying the aim of this workflow [#61](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/61)
- tests on small synthetic genomes + reads [#78](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/78)

### `Fixed`

- multiqc reports of mosdepth and fastp [#63](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/63)
- change fractionBasesAboveQualityThreshold to percentBasesAboveQualityThreshold following meta json file [#63](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/63)
- picard TMP_DIR [#63](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/63)
- reference reading and bug deleting saved bwa reference dir [#73](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/73)

### `Dependencies`

### `Deprecated`

- pytest [#84](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/84)

## v0.2.0-alpha - [19.03.2025]

Initial release of nf-core/grzqc, created with the nf-core template.
