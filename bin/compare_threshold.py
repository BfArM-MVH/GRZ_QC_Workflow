#!/usr/bin/env python3

import pandas as pd
import gzip
import json
import argparse
import sys
from pathlib import Path

def read_bed_file(
    file_path,
    column_names: list | tuple | None = None,
    dtypes: dict[str, str] | None = None,
):
    file_path = Path(file_path)
    is_gzipped = file_path.name.endswith(".gz")

    if column_names is None:
        open_func = gzip.open if is_gzipped else open
        with open_func(file_path, "rb") as f:
            first_line = f.readline().strip().decode("utf-8").split("\t")
            num_columns = len(first_line)

        default_column_names = [
            "chrom", 
            "start", 
            "end", 
            "name", 
            "score", 
            "strand",
            "thickStart", 
            "thickEnd", 
            "itemRgb", 
            "blockCount",
            "blockSizes", 
            "blockStarts",
        ]

        if num_columns > len(default_column_names):
            extra_columns = [
                f"extra_col_{i}" for i in range(num_columns - len(default_column_names))
            ]
            column_names = default_column_names + extra_columns
        else:
            column_names = default_column_names[:num_columns]

    dtype_dict = {
        "chrom": "str", 
        "start": "int64", 
        "end": "int64", 
        "name": "str",
        "score": "float64", 
        "strand": "str", 
        "thickStart": "int64",
        "thickEnd": "int64", 
        "itemRgb": "str", 
        "blockCount": "int64",
        "blockSizes": "str", 
        "blockStarts": "str",
    }
    dtype_dict.update({col: "str" for col in column_names[12:]})
    dtype_dict.update(dtypes or {})

    bed_df = pd.read_csv(
        file_path,
        sep="\t",
        names=column_names,
        dtype=dtype_dict,
        comment="#",
        compression="gzip" if is_gzipped else None,
    )

    return bed_df

def parse_args(args=None):
    Description = "Compare the results with the thresholds."
    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("--mosdepth_global_summary", "-s", required=True)
    parser.add_argument("--mosdepth_target_regions_bed", "-b", required=True)
    parser.add_argument("--thresholds", "-t", required=True)
    parser.add_argument("--fastp_json", "-f", required=True, nargs="+", help="fastp json file(s)")
    parser.add_argument("--sample_id", "-i", required=True)
    parser.add_argument("--labDataName", "-n", required=True)
    parser.add_argument("--libraryType", "-l", required=True)
    parser.add_argument("--sequenceSubtype", "-a", required=True)
    parser.add_argument("--genomicStudySubtype", "-g", required=True)
    parser.add_argument("--output", "-o", required=True)
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    with open(args.thresholds, "r") as f:
        thresholds_data = json.load(f)

    thresholds = None
    for item in thresholds_data:
        if (
            item["libraryType"] == args.libraryType
            and item["sequenceSubtype"] == args.sequenceSubtype
            and item["genomicStudySubtype"] == args.genomicStudySubtype
        ):
            thresholds = item["thresholds"]
            break

    if thresholds is None:
        raise ValueError(
            "No matching thresholds found for meta values: "
            + args.libraryType + ", " + args.sequenceSubtype + ", " + args.genomicStudySubtype
        )

    mean_depth_of_converage_required = float(thresholds["meanDepthOfCoverage"])
    df = pd.read_csv(args.mosdepth_global_summary, sep="\t")
    row_name = "total_region" if args.libraryType in ["panel", "wes", "panel_lr", "wes_lr"] else "total"
    mean_depth_of_coverage = df.loc[df["chrom"] == row_name, "mean"].item()

    quality_threshold = thresholds["percentBasesAboveQualityThreshold"]["qualityThreshold"]
    percent_bases_above_quality_threshold_required = thresholds[
        "percentBasesAboveQualityThreshold"]['percentBasesAbove']

    total_bases = 0
    total_bases_above_quality = 0

    for fastp_json_file in args.fastp_json:
        with open(fastp_json_file, "r") as f:
            fastp_data = json.load(f)
        fastp_filtering_stats = fastp_data["summary"]["before_filtering"]

        if f"q{quality_threshold}_rate" not in fastp_filtering_stats:
            raise ValueError(
                f"'q{quality_threshold}_rate' not found in fastp summary: {fastp_json_file}\n"
                f"-> Could not determine percentBasesAboveQualityThreshold for 'qualityThreshold': {quality_threshold}."
            )

        file_total_bases = fastp_filtering_stats["total_bases"]
        file_q_rate = fastp_filtering_stats[f"q{quality_threshold}_rate"]

        total_bases += file_total_bases
        total_bases_above_quality += file_total_bases * file_q_rate

    if total_bases == 0:
        percent_bases_above_quality_threshold = 0
    else:
        fraction_bases_above_quality_threshold = total_bases_above_quality / total_bases
        percent_bases_above_quality_threshold = fraction_bases_above_quality_threshold * 100

    min_coverage = int(thresholds["targetedRegionsAboveMinCoverage"]["minCoverage"])
    targeted_regions_above_min_coverage_required = float(
        thresholds["targetedRegionsAboveMinCoverage"]["fractionAbove"]
    )

    mosdepth_target_regions_df = read_bed_file(
        args.mosdepth_target_regions_bed,
        column_names=["chrom", "start", "end", "coverage"],
        dtypes={"coverage": "float64"},
    )

    if mosdepth_target_regions_df.empty:
        targeted_regions_above_min_coverage = 0
    else:
        targeted_regions_above_min_coverage = (
            (mosdepth_target_regions_df["coverage"] > min_coverage).mean().item()
        )

    quality_check_passed = (
        mean_depth_of_coverage >= mean_depth_of_converage_required
        and percent_bases_above_quality_threshold >= percent_bases_above_quality_threshold_required
        and targeted_regions_above_min_coverage >= targeted_regions_above_min_coverage_required
    )

    qc_df = pd.DataFrame(
        {
            "sampleId": [args.sample_id],
            "labDataName": [args.labDataName],
            "libraryType": [args.libraryType],
            "sequenceSubtype": [args.sequenceSubtype],
            "genomicStudySubtype": [args.genomicStudySubtype],
            "meanDepthOfCoverage": [mean_depth_of_coverage],
            "meanDepthOfCoverageRequired": [mean_depth_of_converage_required],
            "percentBasesAboveQualityThreshold": [percent_bases_above_quality_threshold],
            "qualityThreshold": [quality_threshold],
            "percentBasesAboveQualityThresholdRequired": [percent_bases_above_quality_threshold_required],
            "targetedRegionsAboveMinCoverage": [targeted_regions_above_min_coverage],
            "minCoverage": [min_coverage],
            "targetedRegionsAboveMinCoverageRequired": [targeted_regions_above_min_coverage_required],
            "passedQC": [quality_check_passed],
        }
    )

    qc_df.to_csv(args.output, index=False)

if __name__ == "__main__":
    sys.exit(main())
