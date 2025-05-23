#!/usr/bin/env python3
import sys
from pathlib import Path
import gzip
import argparse
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
    dtype_dict.update(
        {col: "str" for col in column_names[12:]}
    )  # Default extra columns to string

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

def load_mapping(mapping_file):
    mapping_file = Path(mapping_file)
    # Read the mapping file as a DataFrame
    df = pd.read_csv(mapping_file, sep="\t")
    # Check if the required columns exist
    required_cols = ['UCSC.style.name', 'Sequence.Name']
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Error: Mapping file is missing required column(s): {', '.join(missing)}")
    # Build a mapping from NCBI-style (Sequence.Name) to UCSC-style (UCSC.style.name)
    mapping = dict(zip(df['Sequence.Name'], df['UCSC.style.name']))
    return mapping


def main():
    parser = argparse.ArgumentParser(description="Convert a BED file from NCBI to UCSC style if needed")
    parser.add_argument("--bed", "-b", required=True, help="Input BED file (or gzipped BED file)")
    parser.add_argument("--mapping",  "-m", required=True, help="Mapping file with columns 'UCSC.style.name' and 'Sequence.Name'")
    parser.add_argument("--output", "-o", required=True, help="Output converted BED file")
    
    args = parser.parse_args()
    bed_path = Path(args.bed)
    mapping_path = Path(args.mapping)
    out_path = Path(args.output)

    bed_df = read_bed_file(bed_path)

    # load the mapping and convert the 'chrom' column.
    mapping = load_mapping(mapping_path)
    bed_df["chrom"] = bed_df["chrom"].replace(mapping)

    # Write the (possibly converted) BED file.
    bed_df[["chrom", "start", "end"]].to_csv(out_path, sep="\t", index=False, header=False)
    return 0

if __name__ == "__main__":
    sys.exit(main())
