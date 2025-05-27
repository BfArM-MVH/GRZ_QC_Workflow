#!/usr/bin/env python3

import argparse
import importlib.resources
import json

import grz_pydantic_models.resources
import pandas as pd


def main(args: argparse.Namespace):
    threshold_defs = json.loads(
        (
            importlib.resources.files(grz_pydantic_models.resources) / "thresholds.json"
        ).read_text()
    )
    keys2threshold = {}
    for threshold_def in threshold_defs:
        key = (
            threshold_def["libraryType"],
            threshold_def["sequenceSubtype"],
            threshold_def["genomicStudySubtype"],
        )
        if key in keys2threshold:
            raise ValueError(
                "Thresholds definition file contains duplicate definitions"
            )
        keys2threshold[key] = threshold_def

    thresholds = keys2threshold[
        (args.libraryType, args.sequenceSubtype, args.genomicStudySubtype)
    ]["thresholds"]

    # threshold - mean depth of coverage
    mosdepth_summary_df = pd.read_csv(
        args.mosdepth_global_summary, sep="\t", header=0, index_col="chrom"
    )
    mean_depth_of_coverage = mosdepth_summary_df.loc["total_region", "mean"]
    mean_depth_of_coverage_required = float(thresholds["meanDepthOfCoverage"])

    # percent deviation - mean depth of coverage
    pct_dev_mean_depth_of_coverage = (
        mean_depth_of_coverage - args.meanDepthOfCoverage
    ) / args.meanDepthOfCoverage

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

    # percent deviation - percent bases above quality threshold
    pct_dev_percent_bases_above_quality_threshold = (
        percent_bases_above_quality_threshold - args.percentBasesAboveQualityThreshold
    ) / args.percentBasesAboveQualityThreshold

    # Minimum coverage of target regions to pass
    min_coverage = thresholds["targetedRegionsAboveMinCoverage"]["minCoverage"]
    # Fraction of target regions that must have a coverage above the minimum coverage threshold to pass the validation
    targeted_regions_above_min_coverage_required = thresholds[
        "targetedRegionsAboveMinCoverage"
    ]["fractionAbove"]

    # Read mosdepth target region result
    mosdepth_target_regions_df = pd.read_csv(
        args.mosdepth_target_regions_bed,
        sep="\t",
        names=["chrom", "start", "end", "coverage"],
        usecols=["coverage"],
        dtype={"coverage": float},
    )
    # Compute the fraction of the target regions that have a coverage above the threshold
    if mosdepth_target_regions_df.empty:
        targeted_regions_above_min_coverage = 0
    else:
        targeted_regions_above_min_coverage = (
            mosdepth_target_regions_df["coverage"] >= min_coverage
        ).mean()

    pct_dev_targeted_regions_above_min_coverage = (
        targeted_regions_above_min_coverage - args.targetedRegionsAboveMinCoverage
    ) / args.targetedRegionsAboveMinCoverage

    ### Perform the quality check
    quality_check_passed = (
        (pct_dev_mean_depth_of_coverage <= 0.05)
        and (pct_dev_percent_bases_above_quality_threshold <= 0.05)
        and (pct_dev_targeted_regions_above_min_coverage <= 0.05)
    )

    ### Write the results to a CSV file
    qc_df = pd.DataFrame(
        {
            "sampleId": [args.sample_id],
            "labDataName": [args.labDataName],
            "libraryType": [args.libraryType],
            "sequenceSubtype": [args.sequenceSubtype],
            "genomicStudySubtype": [args.genomicStudySubtype],
            "qualityControlStatus": ["PASS" if quality_check_passed else "FAIL"],
            "meanDepthOfCoverage": [mean_depth_of_coverage],
            "meanDepthOfCoverageRequired": [mean_depth_of_coverage_required],
            "meanDepthOfCoverageDeviation": [pct_dev_mean_depth_of_coverage],
            "percentBasesAboveQualityThreshold": [
                percent_bases_above_quality_threshold
            ],
            "qualityThreshold": [quality_threshold],
            "percentBasesAboveQualityThresholdRequired": [
                percent_bases_above_quality_threshold_required
            ],
            "percentBasesAboveQualityThresholdDeviation": [
                pct_dev_percent_bases_above_quality_threshold
            ],
            "targetedRegionsAboveMinCoverage": [targeted_regions_above_min_coverage],
            "minCoverage": [min_coverage],
            "targetedRegionsAboveMinCoverageRequired": [
                targeted_regions_above_min_coverage_required
            ],
            "targetedRegionsAboveMinCoverageDeviation": [
                pct_dev_targeted_regions_above_min_coverage
            ],
        }
    )
    # write QC results to a CSV file
    qc_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare the results with the thresholds."
    )
    parser.add_argument("--mosdepth_global_summary", "-s", required=True)
    parser.add_argument("--mosdepth_target_regions_bed", "-b", required=True)
    parser.add_argument(
        "--fastp_json", "-f", required=True, nargs="+", help="fastp json file(s)"
    )
    parser.add_argument("--sample_id", "-i", required=True)
    parser.add_argument("--labDataName", "-n", required=True)
    parser.add_argument("--libraryType", "-l", required=True)
    parser.add_argument("--sequenceSubtype", "-a", required=True)
    parser.add_argument("--genomicStudySubtype", "-g", required=True)
    parser.add_argument("--meanDepthOfCoverage", "-m", required=True, type=float)
    parser.add_argument(
        "--targetedRegionsAboveMinCoverage", "-t", required=True, type=float
    )
    parser.add_argument(
        "--percentBasesAboveQualityThreshold", "-p", required=True, type=float
    )
    parser.add_argument("--output", "-o", required=True)
    args = parser.parse_args()

    main(args)
