# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.3.0 - [08.12.2025]

- Derive required thresholds directly from the submission metadata.json / samplesheet.csv [#171](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/171)
  - Use new thresholds API from `grz-pydantic-models` 2.5.0. This will allow to update required thresholds by bumping `grz-pydantic-models` from now on.
  - Add check whether provided metrics pass the required thresholds
- Add symbolic link to final report.csv [#175](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/175)
- Fix conda issues: multiqc is updated to 1.32 as conda was broken in the old version [#173](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/173)
- Disable filters in fastp and fastplong, and report both deduplicated and non-deduplicated reads [#161](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/165)
- Add `micromamba` profile to support `micromamba` as a single-binary conda replacement [#156](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/156)
- Bump grz-pydantic-models to 2.4.0 [#160](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/160)
  - uses new mean read length threshold
- Add instructions for bumping Docker/Singularity local module environments [#162](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/162)

## v1.2.0 - [09.10.2025]

- Ignore non-QCed library types in submission metadata (e.g. RNA) [#149](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/149)
- The `calculate_basequality` module is deprecated as it was not performant for big files. fastp JSON input is not mandatory for aligned file input.
- Bump grz-pydantic-models to support metadata v1.3 [#153](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/153)
- Add parameter to enable sorting paired-end FASTQ files before alignment [#152](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/152)

## v1.1.5 - [01.08.2025]

### Added

- High-level diagram of workflow [#142](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/142)

### Fixed

- Use correct minimap2 prefix for PacBio samples [#147](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/147)

## v1.1.4 - [28.07.2025]

### Fixed

- Fix manifest version to match release tag [#145](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/145)

## v1.1.3 - [25.07.2025]

### Fixed

- Bump grz-pydantic-models to v2.1.2 to allow missing VCFs [#144](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/144)
- Patch bwa-mem2 nf-core modules to use v2.3 to fix possible memory-related crash [#143](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/143)

## v1.1.2 - [17.07.2025]

### Fixed

- Only fail QC if pipeline computes metrics more than 5% less than reported values [#141](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/141)

## v1.1.1 - [15.07.2025]

### Fixed

- Improve QC reports by @twrightsman in [#137](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/137)
- Use latest grz-pydantic-models for end-to-end tests by @twrightsman in [#140](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/140)

## v1.1.0 - [04.07.2025]

### Added

- Support for Nanopore submitted as BAM [#133](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/133)

### Fixed

- Update manifest version to match tag
- Update changelog specification link to latest release

## v1.0.3 - [01.07.2025]

### Fixed

- Fix bug where runs with alignment input are incomplete [#124](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/130)

## v1.0.2 - [23.06.2025]

### Fixed

- Fix bug where fastplong channel is not used as input for COMPARE_THRESHOLDS [#124](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/124)

## v1.0.1 - [20.06.2025]

### Fixed

- Fix bug where only one sample would be aligned due to reference channel being a queue instead of a value [#123](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/123)

## v1.0.0 - [18.06.2025]

### `Added`

- Mark duplicates [#41](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/41)
- decrease reference reading to one (either GRCh37 or GRCh38) [#36](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/36)
- two inputs: metadata path or samplesheet [#54](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/54)
- observed resource usage on README [#35](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/35)
- Save workflow reports from CI [#48](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/48)
- update README about usage of costum config + update base.config for bwamem [#63](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/63)
- Add disclaimer clarifying the aim of this workflow [#61](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/61)
- tests on small synthetic genomes + reads [#78](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/78)
- sambamba [#91](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/91)
- Stating the workflow from alignment is possible [#99](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/99)
- Improve sample naming throughout pipeline [#106](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/106)
- Add support for long read samples to the pipeline [#103](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/103)
- Simplify reference caching to just `--reference_path` [#115](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/115)

### `Fixed`

- multiqc reports of mosdepth and fastp [#63](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/63)
- change fractionBasesAboveQualityThreshold to percentBasesAboveQualityThreshold following meta json file [#63](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/63)
- picard TMP_DIR [#63](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/63)
- reference reading and bug deleting saved bwa reference dir [#73](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/73)
- Migrate GIABv3 reference for GRCh38 [#88](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/88)
- Merging different lanes for the same sample after alignment (BWAMEM2_MEM) instead of merging FASTQs by CAT. Preserving "@RG" is helpfull to track lane based issues. [#90](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/90)
- Use correct QC logic (% deviation) to determine pass/fail [#100](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/100)
- use conda-based package and remove pip for conda envs [#117](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/117)

### `Dependencies`

### `Deprecated`

- pytest [#84](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/84)
- picard [#91](https://github.com/BfArM-MVH/GRZ_QC_Workflow/pull/91)

## v0.2.0-alpha - [19.03.2025]

Initial release of nf-core/grzqc, created with the nf-core template.
