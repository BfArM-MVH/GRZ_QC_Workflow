process CHECK_THRESHOLDS {
  tag "$meta.id"
  label 'process_single'

  input:
    tuple val(meta), path(fastp_json), path(mosdepth_summary_txt), path(mosdepth_region_bed)

  output:
    path('*.threshold_mqc.tsv')

  script:
    """
    #!/usr/bin/env python
    import csv
    import importlib.resources
    import json
    from pathlib import Path
    
    import grz_pydantic_models.resources
    import polars as pl
    
    
    sample_id = "$meta.id"
    lab_data_name = "$meta.lab_data_name"
    library_type = "$meta.library_type"
    sequence_subtype = "$meta.sequence_subtype"
    genomic_study_subtype = "$meta.genomic_study_subtype"
    fastp_stats_path = Path("$fastp_json")
    mosdepth_summary_path = Path("$mosdepth_summary_txt")
    mosdepth_region_bed_path = Path("$mosdepth_region_bed")
    
    threshold_defs = json.loads((importlib.resources.files(grz_pydantic_models.resources) / 'thresholds.json').read_text())
    keys2threshold = {}
    for threshold_def in threshold_defs:
        key = (threshold_def['libraryType'], threshold_def['sequenceSubtype'], threshold_def['genomicStudySubtype'])
        if key in keys2threshold:
            raise ValueError("Thresholds definition file contains duplicate definitions")
        keys2threshold[key] = threshold_def
    
    thresholds = keys2threshold[(library_type, sequence_subtype, genomic_study_subtype)]['thresholds']
    
    mosdepth_summary = pl.read_csv(mosdepth_summary_path, separator = "\t")
    mean_depth_of_coverage = mosdepth_summary.row(
        by_predicate = (pl.col.chrom == 'total_region'),
        named = True
    )['mean']
    pass_coverage_threshold = mean_depth_of_coverage >= thresholds['meanDepthOfCoverage']
    
    mosdepth_region_bed = pl.read_csv(
        source = mosdepth_region_bed_path,
        separator = "\t",
        has_header = False,
        new_columns = ['chrom', 'start', 'end', 'depth']
    )
    fraction_targeted_regions_above_min = mosdepth_region_bed.select(pl.col.depth >= thresholds['targetedRegionsAboveMinCoverage']['minCoverage']).sum().item() / mosdepth_region_bed.height
    pass_target_threshold = fraction_targeted_regions_above_min >= thresholds['targetedRegionsAboveMinCoverage']['fractionAbove']
    
    with open(fastp_stats_path) as fastp_stats_file:
        fastp_stats = json.load(fastp_stats_file)
    quality_threshold = thresholds['percentBasesAboveQualityThreshold']['qualityThreshold']
    pass_quality_threshold = (fastp_stats['summary']['before_filtering'][f"q{quality_threshold}_rate"] * 100) >= thresholds['percentBasesAboveQualityThreshold']['percentBasesAbove']
    
    pass_quality_check = pass_coverage_threshold and pass_target_threshold and pass_quality_threshold
    
    with open(f"{sample_id}.threshold_mqc.tsv", 'w') as output_file:
        output_file.write("\\n".join([
            '# id: "grz_qc"',
            '# section_name: "GRZ QC"',
            '# description: "Results from the GRZ QC threshold check"',
            '# format: "tsv"',
            '# plot_type: "table"'
        ]) + "\\n")
        pl.DataFrame({
            'sample_id': [sample_id],
            'lab_data_name': [lab_data_name],
            'library_type': [library_type],
            'sequence_subtype': [sequence_subtype],
            'genomic_study_subtype': [genomic_study_subtype],
            'quality_check_result': ["PASS" if pass_quality_check else "FAIL"],
            'coverage_depth_mean_observed': [mean_depth_of_coverage],
            'coverage_depth_mean_threshold': [thresholds['meanDepthOfCoverage']],
            'frac_targets_above_min_cov_observed': [fraction_targeted_regions_above_min],
            'frac_targets_above_min_cov_threshold': [thresholds['targetedRegionsAboveMinCoverage']['fractionAbove']],
            'perc_bases_above_min_qual_observed': [fastp_stats['summary']['before_filtering'][f"q{quality_threshold}_rate"] * 100],
            'perc_bases_above_min_qual_threshold': [thresholds['percentBasesAboveQualityThreshold']['percentBasesAbove']],
        }).write_csv(
            file = output_file,
            separator = "\\t"
        )
    """
}
