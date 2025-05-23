#!/usr/bin/env python3

import argparse
import gzip
import json
import sys
from pathlib import Path

import pandas as pd


def read_bed_file(
    file_path,
    column_names: list | tuple | None = None,
    dtypes: dict[str, str] | None = None,
):
    file_path = Path(file_path)
    is_gzipped = file_path.name.endswith(".gz")

    if column_names is None:
        # Determine if the file is compressed
        open_func = gzip.open if is_gzipped else open
        # Load the file and infer the number of columns from the first line
        with open_func(file_path, "rb") as f:
            first_line = f.readline().strip().decode("utf-8").split("\t")
            num_columns = len(first_line)

        # Default column names for the first 12 standard BED fields
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
        # If there are more columns than default names, add generic names
        if num_columns > len(default_column_names):
            extra_columns = [
                f"extra_col_{i}" for i in range(num_columns - len(default_column_names))
            ]
            column_names = default_column_names + extra_columns
        else:
            column_names = default_column_names[:num_columns]

    # Define data types (dtypes) for the columns
    dtype_dict = {
        "chrom": "str",  # Chromosome names are typically strings
        "start": "int64",  # Start position is integer
        "end": "int64",  # End position is integer
        "name": "str",  # Name is typically a string
        "score": "float64",  # Score is usually a float (can also be integer)
        "strand": "str",  # Strand is a string (either '+' or '-')
        "thickStart": "int64",  # thickStart is an integer
        "thickEnd": "int64",  # thickEnd is an integer
        "itemRgb": "str",  # itemRgb is a string (RGB value)
        "blockCount": "int64",  # blockCount is an integer
        "blockSizes": "str",  # blockSizes is a string (comma-separated list)
        "blockStarts": "str",  # blockStarts is a string (comma-separated list)
    }
    # Apply dtypes to extra columns if present
    dtype_dict.update({col: "str" for col in column_names[12:]})
    # Update with user-specified dtypes (overwrites defaults)
    dtype_dict.update(dtypes or {})
    # Read the BED file with inferred column names and dtypes
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
    parser.add_argument(
        "--fastp_json", "-f", required=True, nargs="+", help="fastp json file(s)"
    )
    parser.add_argument("--sample_id", "-i", required=True)
    parser.add_argument("--labDataName", "-n", required=True)
    parser.add_argument("--libraryType", "-l", required=True)
    parser.add_argument("--sequenceSubtype", "-a", required=True)
    parser.add_argument("--genomicStudySubtype", "-g", required=True)
    parser.add_argument("--output", "-o", required=True)
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    # Read thresholds JSON
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
            + args.libraryType
            + ", "
            + args.sequenceSubtype
            + ", "
            + args.genomicStudySubtype
        )

    ### Collect all statistics

    # --- Determine 'meanDepthOfCoverage' ---

    # Required mean depth of coverage to pass the validation
    mean_depth_of_converage_required = float(thresholds["meanDepthOfCoverage"])
    # Read mosdepth summary file
    df = pd.read_csv(args.mosdepth_global_summary, sep="\t")
    row_name = (
        "total_region"
        if args.libraryType in ["panel", "wes", "panel_lr", "wes_lr"]
        else "total"
    )
    mean_depth_of_coverage = df.loc[df["chrom"] == row_name, "mean"].item()

    # --- Determine 'percentBasesAboveQualityThreshold' ---

    # Base quality threshold
    quality_threshold = thresholds["percentBasesAboveQualityThreshold"][
        "qualityThreshold"
    ]
    percent_bases_above_quality_threshold_required = thresholds[
        "percentBasesAboveQualityThreshold"
    ]["percentBasesAbove"]

    total_bases = 0
    total_bases_above_quality = 0

    # read fastp json files
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

        total_bases += file_total_bases
        total_bases_above_quality += fastp_filtering_stats[
            f"q{quality_threshold}_bases"
        ]

    if total_bases == 0:
        percent_bases_above_quality_threshold = 0
    else:
        fraction_bases_above_quality_threshold = total_bases_above_quality / total_bases
        percent_bases_above_quality_threshold = (
            fraction_bases_above_quality_threshold * 100
        )

    # Minimum coverage of target regions to pass
    min_coverage = int(thresholds["targetedRegionsAboveMinCoverage"]["minCoverage"])
    # Fraction of target regions that must have a coverage above the minimum coverage threshold to pass the validation
    targeted_regions_above_min_coverage_required = float(
        thresholds["targetedRegionsAboveMinCoverage"]["fractionAbove"]
    )

    # Read mosdepth target region result
    mosdepth_target_regions_df = read_bed_file(
        args.mosdepth_target_regions_bed,
        column_names=["chrom", "start", "end", "coverage"],
        dtypes={"coverage": "float64"},
    )
    # Compute the fraction of the target regions that have a coverage above the threshold
    if mosdepth_target_regions_df.empty:
        targeted_regions_above_min_coverage = 0
    else:
        targeted_regions_above_min_coverage = (
            (mosdepth_target_regions_df["coverage"] > min_coverage).mean().item()
        )
    ### Perform the quality check
    quality_check_passed = (
        mean_depth_of_coverage >= mean_depth_of_converage_required
        and percent_bases_above_quality_threshold
        >= percent_bases_above_quality_threshold_required
        and targeted_regions_above_min_coverage
        >= targeted_regions_above_min_coverage_required
    )

    ### Write the results to a CSV file
    qc_df = pd.DataFrame(
        {
            "sampleId": [args.sample_id],
            "labDataName": [args.labDataName],
            "libraryType": [args.libraryType],
            "sequenceSubtype": [args.sequenceSubtype],
            "genomicStudySubtype": [args.genomicStudySubtype],
            "meanDepthOfCoverage": [mean_depth_of_coverage],
            "meanDepthOfCoverageRequired": [mean_depth_of_converage_required],
            "percentBasesAboveQualityThreshold": [
                percent_bases_above_quality_threshold
            ],
            "qualityThreshold": [quality_threshold],
            "percentBasesAboveQualityThresholdRequired": [
                percent_bases_above_quality_threshold_required
            ],
            "targetedRegionsAboveMinCoverage": [targeted_regions_above_min_coverage],
            "minCoverage": [min_coverage],
            "targetedRegionsAboveMinCoverageRequired": [
                targeted_regions_above_min_coverage_required
            ],
            "passedQC": [quality_check_passed],
        }
    )
    # write QC results to a CSV file
    qc_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    sys.exit(main())
