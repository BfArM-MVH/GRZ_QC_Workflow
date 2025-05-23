#!/usr/bin/env python3

import pysam
import argparse
import sys
import json


def parse_args(args=None):
    description = "Calculate base quality above the target"

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--output", "-o", required=True, help="Output JSON file path")
    parser.add_argument("--input", "-i", required=True, help="Input alignment file")
   
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    # Open BAM file
    bamfile = pysam.AlignmentFile(args.input, "rb")

    # Initialize counters
    total_reads = 0
    total_bases = 0
    q20_bases = 0
    q30_bases = 0

    for read in bamfile.fetch(until_eof=True):

        quals = read.query_qualities
        if quals is None:
            continue

        total_reads += 1
        total_bases += len(quals)
        q20_bases += sum(q >= 20 for q in quals)
        q30_bases += sum(q >= 30 for q in quals)

    bamfile.close()

    # Avoid division by zero
    q20_rate = q20_bases / total_bases if total_bases > 0 else 0
    q30_rate = q30_bases / total_bases if total_bases > 0 else 0

    # Prepare results
    results = {
        "summary": {
            "before_filtering": {
                "total_reads": total_reads,
                "total_bases": total_bases,
                "q20_bases": q20_bases,
                "q30_bases": q30_bases,
                "q20_rate": round(q20_rate, 6),
                "q30_rate": round(q30_rate, 6)
            }
        }
    }

    # Write to JSON file
    with open(args.output, "w") as f:
        json.dump(results, f, indent=4)

    print(f"Results written to {args.output}")

if __name__ == "__main__":
    sys.exit(main())
