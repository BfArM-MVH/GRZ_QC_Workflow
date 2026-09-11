#!/usr/bin/env python3

import argparse
import json
from enum import StrEnum

import pandas as pd

PCT_DEV_CUTOFF = 10


class QCStatus(StrEnum):
    PASS = "PASS"
    DEVIATION = "DEVIATION"
    THRESHOLD_NOT_MET = "THRESHOLD NOT MET"


def percent_deviation(computed: float, provided: float | None) -> float | None:
    """Signed percent deviation between the computed and provided values, relative to the
    computed value (the one determined by the GRZ).

    The result is negative when the provided value exceeds the computed value (i.e. the
    Leistungserbringer reported too high) and positive when it is below.

    :param computed: value computed by this pipeline
    :param provided: value provided in the submission metadata, if any
    :return: signed percent deviation, or ``None`` if no value was provided or the computed
        value is zero (deviation undefined)
    """
    # A provided value of zero is treated like a missing one (the Nextflow module already
    # omits zero-valued provided metrics); a computed value of zero leaves the relative
    # deviation undefined.
    if provided is None or provided == 0 or computed == 0:
        return None
    return (computed - provided) / computed * 100


def classify_metric(
    threshold_passed: bool, deviation_percent: float | None
) -> QCStatus:
    """Determine the QC status of a single metric.

    :param threshold_passed: whether the value computed by the pipeline meets the required
        BfArM threshold
    :param deviation_percent: percent deviation between the computed and provided
        values (see :func:`percent_deviation`), or ``None`` if no value was provided
    :return: ``THRESHOLD NOT MET`` if the computed value is below the required threshold
        (fails QC); ``DEVIATION`` if the provided value deviates by more than
        ``PCT_DEV_CUTOFF`` percent from the computed value in either direction (reported, but
        does not fail QC); ``PASS`` otherwise
    """
    if not threshold_passed:
        return QCStatus.THRESHOLD_NOT_MET
    if deviation_percent is not None and abs(deviation_percent) > PCT_DEV_CUTOFF:
        return QCStatus.DEVIATION
    return QCStatus.PASS


def main(args: argparse.Namespace):
    # threshold - mean depth of coverage
    mosdepth_summary_df = pd.read_csv(
        args.mosdepth_global_summary, sep="\t", header=0, index_col="chrom"
    )

    if args.libraryType.startswith("wgs"):
        # WGS should compute mean depth over entire genome
        depth_key = "total"
    elif args.libraryType.split("_")[0] in ["panel", "wes"]:
        # panel + WES only over the specified regions
        depth_key = "total_region"
    else:
        raise ValueError(f"Unknown library type: {args.libraryType}")

    mean_depth_of_coverage_provided = args.meanDepthOfCoverage
    mean_depth_of_coverage_computed = mosdepth_summary_df.loc[depth_key, "mean"]
    mean_depth_of_coverage_required = float(args.meanDepthOfCoverageRequired)

    # threshold check - mean depth of coverage
    mean_depth_of_coverage__threshold_passed = (
        mean_depth_of_coverage_computed >= mean_depth_of_coverage_required
    )

    # percent deviation - mean depth of coverage
    mean_depth_of_coverage__deviation_percent = percent_deviation(
        mean_depth_of_coverage_computed, mean_depth_of_coverage_provided
    )

    # Base quality threshold
    quality_threshold = args.qualityThreshold
    percent_bases_above_quality_threshold_provided = (
        args.percentBasesAboveQualityThreshold
    )
    percent_bases_above_quality_threshold_required = (
        args.percentBasesAboveQualityThresholdRequired
    )

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
        percent_bases_above_quality_threshold_computed = 0
    else:
        fraction_bases_above_quality_threshold_computed = (
            total_bases_above_quality / total_bases
        )
        percent_bases_above_quality_threshold_computed = (
            fraction_bases_above_quality_threshold_computed * 100
        )

    # threshold check - percent bases above quality threshold
    percent_bases_above_quality_threshold__threshold_passed = (
        percent_bases_above_quality_threshold_computed
        >= percent_bases_above_quality_threshold_required
    )

    # percent deviation - percent bases above quality threshold
    percent_bases_above_quality_threshold__deviation_percent = percent_deviation(
        percent_bases_above_quality_threshold_computed,
        percent_bases_above_quality_threshold_provided,
    )

    # Minimum coverage of target regions to pass
    min_coverage = args.minCoverage
    targeted_regions_above_min_coverage_provided = args.targetedRegionsAboveMinCoverage
    # Fraction of target regions that must have a coverage above the minimum coverage threshold to pass the validation
    targeted_regions_above_min_coverage_required = (
        args.targetedRegionsAboveMinCoverageRequired
    )

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
        targeted_regions_above_min_coverage_computed = 0
    else:
        targeted_regions_above_min_coverage_computed = (
            mosdepth_target_regions_df["coverage"] >= min_coverage
        ).mean()

    # threshold check - targeted regions above minimum coverage
    targeted_regions_above_min_coverage__threshold_passed = (
        targeted_regions_above_min_coverage_computed
        >= targeted_regions_above_min_coverage_required
    )

    # percent deviation - targeted regions above minimum coverage
    targeted_regions_above_min_coverage__deviation_percent = percent_deviation(
        targeted_regions_above_min_coverage_computed,
        targeted_regions_above_min_coverage_provided,
    )

    ### Perform the quality check(s)
    mean_depth_of_coverage__qc_status = classify_metric(
        mean_depth_of_coverage__threshold_passed,
        mean_depth_of_coverage__deviation_percent,
    )
    percent_bases_above_quality_threshold__qc_status = classify_metric(
        percent_bases_above_quality_threshold__threshold_passed,
        percent_bases_above_quality_threshold__deviation_percent,
    )
    targeted_regions_above_min_coverage__qc_status = classify_metric(
        targeted_regions_above_min_coverage__threshold_passed,
        targeted_regions_above_min_coverage__deviation_percent,
    )

    # A metric flagged DEVIATION is reported to the Plattformträger but does not fail the
    # Detailprüfung, as long as the value computed here still meets the BfArM threshold.
    quality_check_passed = QCStatus.THRESHOLD_NOT_MET not in (
        mean_depth_of_coverage__qc_status,
        percent_bases_above_quality_threshold__qc_status,
        targeted_regions_above_min_coverage__qc_status,
    )
    quality_control_status = "PASS" if quality_check_passed else "FAIL"

    ### Write the results to a CSV file
    qc_df = pd.DataFrame(
        {
            "sampleId": [args.sample_id],
            "donorPseudonym": [args.donorPseudonym],
            "labDataName": [args.labDataName],
            "libraryType": [args.libraryType],
            "sequenceSubtype": [args.sequenceSubtype],
            "genomicStudySubtype": [args.genomicStudySubtype],
            "qualityControlStatus": [quality_control_status],
            "meanDepthOfCoverage": [mean_depth_of_coverage_computed],
            "meanDepthOfCoverageProvided": [args.meanDepthOfCoverage],
            "meanDepthOfCoverageRequired": [mean_depth_of_coverage_required],
            "meanDepthOfCoverageDeviation": [mean_depth_of_coverage__deviation_percent],
            "meanDepthOfCoverageQCStatus": [mean_depth_of_coverage__qc_status],
            "percentBasesAboveQualityThreshold": [
                percent_bases_above_quality_threshold_computed
            ],
            "qualityThreshold": [quality_threshold],
            "percentBasesAboveQualityThresholdProvided": [
                args.percentBasesAboveQualityThreshold
            ],
            "percentBasesAboveQualityThresholdRequired": [
                percent_bases_above_quality_threshold_required
            ],
            "percentBasesAboveQualityThresholdDeviation": [
                percent_bases_above_quality_threshold__deviation_percent
            ],
            "percentBasesAboveQualityThresholdQCStatus": [
                percent_bases_above_quality_threshold__qc_status
            ],
            "targetedRegionsAboveMinCoverage": [
                targeted_regions_above_min_coverage_computed
            ],
            "minCoverage": [min_coverage],
            "targetedRegionsAboveMinCoverageProvided": [
                args.targetedRegionsAboveMinCoverage
            ],
            "targetedRegionsAboveMinCoverageRequired": [
                targeted_regions_above_min_coverage_required
            ],
            "targetedRegionsAboveMinCoverageDeviation": [
                targeted_regions_above_min_coverage__deviation_percent
            ],
            "targetedRegionsAboveMinCoverageQCStatus": [
                targeted_regions_above_min_coverage__qc_status
            ],
        }
    )
    # write QC results to a CSV file
    qc_df.to_csv(args.output, index=False)


def non_negative_float(value: str) -> float:
    """Argparse type: float >= 0, else error matching schema message."""
    try:
        f = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a positive number.")
    if f < 0:
        raise argparse.ArgumentTypeError("Must be a positive number.")
    return f


def fraction_0_1(value: str) -> float:
    """Argparse type: 0 <= float <= 1, else error matching schema message."""
    try:
        f = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a number between 0 and 1.")
    if f < 0 or f > 1:
        raise argparse.ArgumentTypeError("Must be a number between 0 and 1.")
    return f


def percentage_0_100(value: str) -> float:
    """Argparse type: 0 <= float <= 100, else error matching schema message."""
    try:
        f = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a number between 0 and 100.")
    if f < 0 or f > 100:
        raise argparse.ArgumentTypeError("Must be a number between 0 and 100.")
    return f


def non_negative_int(value: str) -> int:
    """Argparse type: int >= 0, else error matching schema message."""
    try:
        i = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a positive integer.")
    if i < 0:
        raise argparse.ArgumentTypeError("Must be a positive integer.")
    return i


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
    parser.add_argument("--donorPseudonym", "-d", required=True)
    parser.add_argument("--libraryType", "-l", required=True)
    parser.add_argument("--sequenceSubtype", "-a", required=True)
    parser.add_argument("--genomicStudySubtype", "-g", required=True)
    parser.add_argument("--meanDepthOfCoverage", "-m", type=non_negative_float)
    parser.add_argument("--targetedRegionsAboveMinCoverage", "-t", type=fraction_0_1)
    parser.add_argument(
        "--percentBasesAboveQualityThreshold", "-p", type=percentage_0_100
    )

    # required thresholds
    parser.add_argument(
        "--meanDepthOfCoverageRequired",
        type=non_negative_float,
        required=True,
        help="Required mean depth of coverage",
    )
    parser.add_argument(
        "--qualityThreshold",
        type=non_negative_int,
        required=True,
        help="Base quality threshold (Q score) used to compute percent of bases above this quality",
    )
    parser.add_argument(
        "--percentBasesAboveQualityThresholdRequired",
        type=percentage_0_100,
        required=True,
        help="Required percent of bases above the quality threshold",
    )
    parser.add_argument(
        "--minCoverage",
        type=non_negative_float,
        required=True,
        help="Minimum per-region coverage threshold for target regions",
    )
    parser.add_argument(
        "--targetedRegionsAboveMinCoverageRequired",
        type=fraction_0_1,
        required=True,
        help="Required fraction of target regions that must be above minCoverage to pass QC",
    )

    parser.add_argument("--output", "-o", required=True)
    args = parser.parse_args()

    main(args)
