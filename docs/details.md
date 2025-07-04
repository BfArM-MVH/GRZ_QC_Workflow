# Pipeline Details

## Mosdepth

[Mosdepth](https://github.com/brentp/mosdepth) reports information for the evaluation of the quality of the provided alignment data.
In short, the basic statistics of the alignment (number of reads, coverage, GC-content, etc.) are summarized and a number of useful graphs are produced.

The pipeline uses branching logic to provide the correct BED file depending on the library type:

- For WGS samples:
  The pipeline automatically selects the default BED file containing ~400 genes for the mosdepth run. Mosdepth calculates depths for both genome-wide and targeted regions in one run.

- For WES/Panel samples:
  The pipeline uses the target region BED file supplied in the submission metadata and converts it into UCSC-style if the chromosome names are in NCBI-style.

| **Library Type** | **input: BED File**                                                                                                                              | `meanDepthOfCoverage`                      | `targetedRegionsAboveMinCoverage`                   | `percentBasesAboveQualityThreshold`                                                          |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------ | --------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| **WGS**          | Default representative genes, ([hg19](../assets/default_files/hg19_439_omim_genes.bed), [hg38](../assets/default_files/hg38_440_omim_genes.bed)) | genome-wide average coverage               | Computed over a predefined set of ~400 gene regions | Proportion of bases in untrimmed+unfiltered reads above minimum quality for the library type |
| **WES/Panel**    | User-provided BED file                                                                                                                           | Average coverage on user-specified targets | Computed over the user-provided target regions      | Proportion of bases in untrimmed+unfiltered reads above minimum quality for the library type |

## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.
