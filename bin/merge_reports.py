#!/usr/bin/env python3

import argparse
from enum import StrEnum
from textwrap import dedent

import pandas as pd

# Keep the following in sync with compare_threshold.py
PCT_DEV_CUTOFF = 10


class QCStatus(StrEnum):
    PASS = "PASS"
    DEVIATION = "DEVIATION"
    THRESHOLD_NOT_MET = "THRESHOLD NOT MET"


def qc_status_description(metric: str) -> str:
    """Human-readable description of a per-metric QC status column."""
    return (
        f"'{QCStatus.THRESHOLD_NOT_MET}' if the computed {metric} is below the value "
        f"required by BfArM (fails QC), '{QCStatus.DEVIATION}' if the provided value "
        f"deviates by more than {PCT_DEV_CUTOFF}% from it in either direction (reported, "
        f"but does not fail QC), '{QCStatus.PASS}' otherwise."
    )


def main(args: argparse.Namespace):
    # concat all csv files
    dfs = [pd.read_csv(f) for f in args.inputs]
    df_merged = pd.concat(dfs, ignore_index=True)
    df_merged["grzQcWorkflowVersion"] = args.workflow_version
    df_merged.to_csv(f"{args.output_prefix}.csv", index=False)
    df_merged.to_excel(f"{args.output_prefix}.xlsx", index=False)

    # write out annotated report for MultiQC
    mean_depth_qc_desc = qc_status_description("mean depth of coverage")
    percent_bases_qc_desc = qc_status_description(
        "percentage of bases above the quality threshold"
    )
    targeted_regions_qc_desc = qc_status_description(
        "proportion of target regions above the minimum coverage"
    )
    with open(f"{args.output_prefix}_mqc.csv", "w") as mqc_out:
        mqc_out.write(
            dedent(f"""\
        # id: "{args.mqc_id}"
        # section_name: "{args.mqc_section_name}"
        # description: "Results from the GRZ internal QC pipeline."
        # format: "csv"
        # plot_type: "table"
        # headers:
        #   donorPseudonym:
        #     title: "Donor Pseudonym"
        #     description: "A unique identifier given by the Leistungserbringer for each donor."
        #   labDataName:
        #     title: "Lab Data Name"
        #     description: "Name/ID of the biospecimen."
        #   libraryType:
        #     title: "Library Type"
        #     description: "Type of sequencing library (e.g. 'wgs')."
        #   sequenceSubtype:
        #     title: "Sequence Subtype"
        #     description: "Subtype of sequence (germline, somatic, etc.)."
        #   genomicStudySubtype:
        #     title: "Genomic Study Subtype"
        #     description: "Whether tumor and/or germline are tested."
        #   qualityControlStatus:
        #     title: "Overall QC Status"
        #     description: "Whether all pipeline-computed metrics meet the thresholds required by BfArM. A provided metric deviating more than {PCT_DEV_CUTOFF}% from the computed value leads to the reporting of the submission but does not lead to a failed quality control."
        #     cond_formatting_rules:
        #       pass:
        #         - s_eq: "PASS"
        #       fail:
        #         - s_eq: "FAIL"
        #   meanDepthOfCoverage:
        #     title: "Mean Depth of Coverage"
        #     description: "Mean depth of coverage computed by the pipeline."
        #   meanDepthOfCoverageProvided:
        #     title: "Mean Depth of Coverage Provided"
        #     description: "Mean depth of coverage provided in the samplesheet/submission metadata."
        #   meanDepthOfCoverageRequired:
        #     title: "Mean Depth of Coverage Required"
        #     description: "Mean depth of coverage required to pass quality control."
        #   meanDepthOfCoverageDeviation:
        #     title: "Mean Depth of Coverage Deviation"
        #     description: "Signed percent deviation between the computed and provided mean depth of coverage, relative to the computed value. Negative if the provided value is higher (the Leistungserbringer reported too high)."
        #     suffix: '%'
        #   meanDepthOfCoverageQCStatus:
        #     title: "Mean Depth of Coverage QC Status"
        #     description: "{mean_depth_qc_desc}"
        #     cond_formatting_rules:
        #       pass:
        #         - s_eq: "{QCStatus.PASS}"
        #       warn:
        #         - s_eq: "{QCStatus.DEVIATION}"
        #       fail:
        #         - s_eq: "{QCStatus.THRESHOLD_NOT_MET}"
        #   percentBasesAboveQualityThreshold:
        #     title: "Percent Bases Above Quality Threshold"
        #     description: "Percentage of unfiltered read bases that are above the minimum quality score."
        #     suffix: '%'
        #   qualityThreshold:
        #     title: "Quality Threshold"
        #     description: "Minimum quality score for computing percentage of bases above it."
        #   percentBasesAboveQualityThresholdProvided:
        #     title: "Percent Bases Above Quality Threshold Provided"
        #     description: "Percentage of unfiltered read bases above minimum quality score provided in the samplesheet/submission metadata."
        #     suffix: '%'
        #   percentBasesAboveQualityThresholdRequired:
        #     title: "Percent Bases Above Quality Threshold Required"
        #     description: "Minimum percentage of unfiltered read bases above minimum quality score to pass quality control."
        #     suffix: '%'
        #   percentBasesAboveQualityThresholdDeviation:
        #     title: "Percent Bases Above Quality Threshold Deviation"
        #     description: "Signed percent deviation between the computed and provided percentage of bases above quality threshold, relative to the computed value. Negative if the provided value is higher (the Leistungserbringer reported too high)."
        #     suffix: '%'
        #   percentBasesAboveQualityThresholdQCStatus:
        #     title: "Percent Bases Above Quality Threshold QC Status"
        #     description: "{percent_bases_qc_desc}"
        #     cond_formatting_rules:
        #       pass:
        #         - s_eq: "{QCStatus.PASS}"
        #       warn:
        #         - s_eq: "{QCStatus.DEVIATION}"
        #       fail:
        #         - s_eq: "{QCStatus.THRESHOLD_NOT_MET}"
        #   targetedRegionsAboveMinCoverage:
        #     title: "Targeted Regions Above Minimum Coverage"
        #     description: "Proportion of target regions above the minimum coverage threshold."
        #   minCoverage:
        #     title: "Targeted Region Minimum Coverage"
        #     description: "Minimum coverage for computing proportion of targeted regions above it."
        #   targetedRegionsAboveMinCoverageProvided:
        #     title: "Targeted Regions Above Minimum Coverage Provided"
        #     description: "Proportion of target regions above the minimum coverage threshold provided in the samplesheet/submission metadata."
        #   targetedRegionsAboveMinCoverageRequired:
        #     title: "Targeted Regions Above Minimum Coverage Required"
        #     description: "Minimum proportion of target regions above the minimum coverage threshold required to pass quality control."
        #   targetedRegionsAboveMinCoverageDeviation:
        #     title: "Targeted Regions Above Minimum Coverage Deviation"
        #     description: "Signed percent deviation between the computed and provided proportion of target regions above the minimum coverage threshold, relative to the computed value. Negative if the provided value is higher (the Leistungserbringer reported too high)."
        #     suffix: '%'
        #   targetedRegionsAboveMinCoverageQCStatus:
        #     title: "Targeted Regions Above Minimum Coverage QC Status"
        #     description: "{targeted_regions_qc_desc}"
        #     cond_formatting_rules:
        #       pass:
        #         - s_eq: "{QCStatus.PASS}"
        #       warn:
        #         - s_eq: "{QCStatus.DEVIATION}"
        #       fail:
        #         - s_eq: "{QCStatus.THRESHOLD_NOT_MET}"
        #   grzQcWorkflowVersion:
        #     title: "GRZ QC Workflow Version"
        #     description: "Version of the GRZ QC workflow used to produce this report."
        """)
        )
        df_merged.to_csv(mqc_out, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare the results with the thresholds."
    )
    parser.add_argument("inputs", nargs="+", help="List of files to merge")
    parser.add_argument(
        "--output_prefix", "-o", required=True, help="Output file prefix"
    )
    parser.add_argument(
        "--workflow_version",
        required=True,
        help="GRZ QC workflow version to record in the report",
    )
    parser.add_argument(
        "--mqc_id",
        default="grz_qc",
        help="MultiQC custom-content section id (must be unique per table to avoid collisions)",
    )
    parser.add_argument(
        "--mqc_section_name",
        default="GRZ QC Results",
        help="MultiQC custom-content section display name",
    )
    args = parser.parse_args()

    main(args)
