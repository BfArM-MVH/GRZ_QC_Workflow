#!/usr/bin/env python3

import argparse
import sys

import pandas as pd


def parse_args(args=None):
    Description = "Compare the results with the thresholds."

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("inputs", nargs="+", help="List of files to merge")
    parser.add_argument(
        "--output_prefix", "-o", required=True, help="Output file prefix"
    )

    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    # concat all csv files
    dfs = [pd.read_csv(f) for f in args.inputs]
    df_merged = pd.concat(dfs, ignore_index=True)
    df_merged.to_csv(f"{args.output_prefix}.csv", index=False)
    df_merged.to_excel(f"{args.output_prefix}.xlsx", index=False)


if __name__ == "__main__":
    sys.exit(main())
