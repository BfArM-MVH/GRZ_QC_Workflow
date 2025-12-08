#!/usr/bin/env python3

import argparse
from textwrap import dedent

import pandas as pd


def main(args: argparse.Namespace):
    # concat all csv files
    dfs = [pd.read_csv(f) for f in args.inputs]
    df_merged = pd.concat(dfs, ignore_index=True)
    df_merged.to_csv(f"{args.output_prefix}.csv", index=False)
    df_merged.to_excel(f"{args.output_prefix}.xlsx", index=False)

    # write out annotated report for MultiQC
    with open(f"{args.output_prefix}_mqc.csv", "w") as mqc_out:
        mqc_out.write(
            dedent("""\
        # id: "grz_qc"
        # section_name: "GRZ QC Results"
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
        #     description: "If pre-computed metrics were provided, this states whether deviation of pipeline-computed metrics from them are all less than the official deviation threshold."
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
        #     description: "Percent deviation of pipeline-computed mean depth of coverage from provided value."
        #     suffix: '%'
        #   meanDepthOfCoverageQCStatus:
        #     title: "Mean Depth of Coverage QC Status"
        #     description: "Whether the sample passes the mean depth of coverage QC criteria or is too low/too high."
        #     cond_formatting_rules:
        #       pass:
        #         - s_eq: "PASS"
        #       fail:
        #         - s_eq: "TOO LOW"
        #         - s_eq: "THRESHOLD NOT MET"
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
        #     description: "Percent deviation of pipeline-computed percentage of bases above quality threshold from provided value."
        #     suffix: '%'
        #   percentBasesAboveQualityThresholdQCStatus:
        #     title: "Percent Bases Above Quality Threshold QC Status"
        #     description: "Whether the sample passes the percent bases above quality threshold QC criteria or is too low/too high."
        #     cond_formatting_rules:
        #       pass:
        #         - s_eq: "PASS"
        #       fail:
        #         - s_eq: "TOO LOW"
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
        #     description: "Percent deviation of pipeline-computed proportion of target regions above the minimum coverage threshold from provided value."
        #     suffix: '%'
        #   targetedRegionsAboveMinCoverageQCStatus:
        #     title: "Targeted Regions Above Minimum Coverage QC Status"
        #     description: "Whether the sample passes the targeted regions above minimum coverage QC criteria or is too low/too high."
        #     cond_formatting_rules:
        #       pass:
        #         - s_eq: "PASS"
        #       fail:
        #         - s_eq: "TOO LOW"
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
    args = parser.parse_args()

    main(args)
