#!/usr/bin/env python3
import pandas as pd
import gzip
import json
import argparse
import sys

def parse_args(args=None):

    Description = "Compare the results with the thresholds."

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("--summary", "-s", required=True, help="mosdepth summary file")
    parser.add_argument("--bed", "-b", required=True, help="mosdepth target bed file")
    parser.add_argument("--thresholds", "-t", required=True, help="thresholds json file")
    parser.add_argument("--fastp_json", "-f", required=True, help="fastp json file")
    parser.add_argument("--sample_id", "-i", required=True, help="sample id/meta.id")
    parser.add_argument("--libraryType", "-l", required=True, help="libraryType")
    parser.add_argument("--sequenceSubtype", "-a", required=True, help="sequenceSubtype")
    parser.add_argument("--genomicStudySubtype", "-g", required=True, help="genomicStudySubtype")

    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    ### Extract results and thresholds
    # Extract Q30 rate from the fastp JSON.
    with open(args.fastp_json, "r") as f:
        fastp_data = json.load(f)
    q30_rate = fastp_data["summary"]["before_filtering"]["q30_rate"]

    # Extract thresholds JSON.
    with open(args.thresholds, "r") as f:
        thresholds_data = json.load(f)

    THRESHOLDS = None
    for item in thresholds_data:
        if (item["libraryType"] == args.libraryType and
            item["sequenceSubtype"] == args.sequenceSubtype and
            item["genomicStudySubtype"] == args.genomicStudySubtype):
            THRESHOLDS = item["thresholds"]
            break

    if THRESHOLDS is None:
        raise ValueError("No matching thresholds found for meta values: " + args.libraryType + ", " + args.sequenceSubtype + ", " + args.genomicStudySubtype)

    MEAN_DEPTH_THRESHOLD = float(THRESHOLDS["meanDepthOfCoverage"])
    QUALITY_THRESHOLD = THRESHOLDS["fractionBasesAboveQualityThreshold"]["qualityThreshold"]

    ## what is default Q30_THRESHOLD value?
    Q30_THRESHOLD = 0.0
    if str(QUALITY_THRESHOLD) == "30":
        Q30_THRESHOLD = float(THRESHOLDS["fractionBasesAboveQualityThreshold"]["fractionBasesAbove"])

    TARGET_MIN_COVERAGE = float(THRESHOLDS["targetedRegionsAboveMinCoverage"]["minCoverage"])
    TARGET_FRACTION_ABOVE_THRESHOLD = float(THRESHOLDS["targetedRegionsAboveMinCoverage"]["fractionAbove"])

    ### compare results with thresholds

    # 1. Read mosdepth summary file
    df = pd.read_csv(args.summary, sep="\\t")
    row_name = "total_region" if args.libraryType in ["panel", "wes"] else "total"
    mosdepth_cov = float(df.loc[df.iloc[:,0] == row_name, "mean"].values[0])

    # 2. parse the mosdepth target gene result
    count = 0
    total = 0
    with gzip.open(args.bed, "rt") as f:
        for line in f:
            total += 1
            parts = line.strip().split()
            if float(parts[3]) >= TARGET_MIN_COVERAGE:
                count += 1

    mosdepth_cov_rate_target = count / total if total > 0 else 0
    ##    # per base solution
    ##    # ! note if we change this , we also need to change the output channel of MOSDEPTH_TARGET
    ##    total_length = 0
    ##    filtered_length = 0
    ##    with gzip.open(bed_file, "rt") as f:
    ##        for line in f:
    ##            parts = line.strip().split()
    ##            # Calculate interval length: col3 - col2
    ##            start = float(parts[1])
    ##            end   = float(parts[2])
    ##            interval_length = end - start
    ##            total_length += interval_length
    ##            # If column 4 > TARGET_MIN_COVERAGE, add to filtered length.
    ##            if float(parts[3]) > TARGET_MIN_COVERAGE:
    ##                filtered_length += interval_length
    ##
    ##    mosdepth_cov_rate_target = filtered_length / total_length if total_length > 0 else 0

    # 3. Perform the quality check.
    quality_check = "PASS" if (mosdepth_cov >= MEAN_DEPTH_THRESHOLD and q30_rate >= Q30_THRESHOLD and mosdepth_cov_rate_target >= TARGET_FRACTION_ABOVE_THRESHOLD) else "FAIL"

    # 4. Write the results to a CSV file.
    with open(f"{args.sample_id}.result.csv", "w") as f:
        f.write("Sample_id,libraryType,sequenceSubtype,genomicStudySubtype,q30_rate,Q30_THRESHOLD,Mosdepth_cov,MEAN_DEPTH_THRESHOLD,Mosdepth_cov_ratio_target_genes,TARGET_FRACTION_ABOVE_THRESHOLD,Quality_check\n")
        f.write(f"{args.sample_id},{args.libraryType},{args.sequenceSubtype},{args.genomicStudySubtype},{q30_rate},{Q30_THRESHOLD},{mosdepth_cov},{MEAN_DEPTH_THRESHOLD},{mosdepth_cov_rate_target},{TARGET_FRACTION_ABOVE_THRESHOLD},{quality_check}\n")

if __name__ == "__main__":
	sys.exit(main())
