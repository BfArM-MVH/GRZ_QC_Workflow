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
