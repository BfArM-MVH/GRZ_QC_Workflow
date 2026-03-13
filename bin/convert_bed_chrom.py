#!/usr/bin/env python3
import argparse
import gzip
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
    # Determine if the file is compressed
    open_func = gzip.open if is_gzipped else open

    # Load the file, identify where data starts and infer the number of columns
    data_start_line = 0
    num_columns = 0
    with open_func(file_path, "rb") as f:
        for line in f:
            line_decoded = line.strip().decode("utf-8")
            if not line_decoded or any(
                line_decoded.startswith(x) for x in ["#", "track", "browser"]
            ):
                data_start_line += 1
                continue
            num_columns = len(line_decoded.split("\t"))
            break

    # Enforce a minimum of three columns
    if num_columns < 3:
        raise ValueError(
            f"Invalid BED: {file_path} has {num_columns} columns. Expected >= 3."
        )

    if column_names is None:
        # Default column names for the first 6 standard BED fields
        default_names = ["chrom", "start", "end", "name", "score", "strand"]

        if num_columns <= len(default_names):
            column_names = default_names[:num_columns]
        else:
            # Add generic names for columns beyond the standard 6
            extra = [f"extra_col_{i}" for i in range(num_columns - len(default_names))]
            column_names = default_names + extra

    # Create dtype_dict dynamically based on columns that exist
    # Use 'Int64' for start/end to allow the strict null check later on
    base_types = {
        "chrom": "str",
        "start": "Int64",
        "end": "Int64",
        "name": "str",
        "score": "float64",
        "strand": "str",
    }

    # Only include types for columns actually present in column_names
    dtype_dict = {col: base_types.get(col, "str") for col in column_names}
    dtype_dict.update(dtypes or {})

    bed_df = pd.read_csv(
        file_path,
        sep="\t",
        names=column_names,
        header=None,
        skiprows=data_start_line,
        dtype=dtype_dict,
        comment="#",
        na_values=["."],
        compression="gzip" if is_gzipped else None,
    )

    # Null Check: Fail if '.' present in chrom, start, or end
    for col in ["chrom", "start", "end"]:
        if bed_df[col].isnull().any():
            row_idx = bed_df[bed_df[col].isnull()].index[0] + data_start_line + 1
            raise ValueError(
                f"ERROR: Null value detected in '{col}' at line {row_idx} of {file_path.name}."
            )

    # Cast to standard int for output
    bed_df["start"] = bed_df["start"].astype(int)
    bed_df["end"] = bed_df["end"].astype(int)

    return bed_df


def load_mapping(mapping_file):
    mapping_file = Path(mapping_file)
    # Read the mapping file as a DataFrame
    df = pd.read_csv(mapping_file, sep="\t")
    # Check if the required columns exist
    required_cols = ["UCSC.style.name", "Sequence.Name"]
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(
            f"Error: Mapping file is missing required column(s): {', '.join(missing)}"
        )
    # Build a mapping from NCBI-style (Sequence.Name) to UCSC-style (UCSC.style.name)
    mapping = dict(zip(df["Sequence.Name"], df["UCSC.style.name"]))
    return mapping


def main():
    parser = argparse.ArgumentParser(
        description="Convert a BED file from NCBI to UCSC style if needed"
    )
    parser.add_argument(
        "--bed", "-b", required=True, help="Input BED file (or gzipped BED file)"
    )
    parser.add_argument(
        "--mapping",
        "-m",
        required=True,
        help="Mapping file with columns 'UCSC.style.name' and 'Sequence.Name'",
    )
    parser.add_argument(
        "--output", "-o", required=True, help="Output converted BED file"
    )

    args = parser.parse_args()
    bed_path = Path(args.bed)
    mapping_path = Path(args.mapping)
    out_path = Path(args.output)

    bed_df = read_bed_file(bed_path)

    # load the mapping and convert the 'chrom' column.
    mapping = load_mapping(mapping_path)
    bed_df["chrom"] = bed_df["chrom"].replace(mapping)

    # Write the (possibly converted) BED file.
    bed_df[["chrom", "start", "end"]].to_csv(
        out_path, sep="\t", index=False, header=False
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
