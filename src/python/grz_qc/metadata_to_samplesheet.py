#!/usr/bin/env python3


import argparse
import pandas as pd
import json
from itertools import groupby
from typing import Dict, Any, Generator

from pathlib import Path


def determine_fastq_pairs(fastq_files: list[dict]) -> list[tuple[dict, dict]]:
    key = lambda f: (f.get("flowcellId", ""), f.get("laneId", ""))

    retval = []

    fastq_files.sort(key=key)
    for _key, group in groupby(fastq_files, key):
        files = list(group)

        # separate R1 and R2 files
        fastq_r1_files = [f for f in files if f.get("readOrder") == "R1"]
        fastq_r2_files = [f for f in files if f.get("readOrder") == "R2"]

        assert len(fastq_r1_files) == 1, (
            f"Expected one R1 file, but got {len(fastq_r1_files)}"
        )
        assert len(fastq_r2_files) == 1, (
            f"Expected one R2 file, but got {len(fastq_r2_files)}"
        )

        retval.append((fastq_r1_files[0], fastq_r2_files[0]))

    return retval


def extract_data(
    json_data: Dict, submission_base_path: str | Path
) -> Generator[dict[str, Path | str | Any], Any, list[Any]]:
    """
    Extract specific fields from the GRZ schema JSON data.

    Args:
        json_data (Dict): The parsed JSON data from the GRZ schema file.
        submission_base_path : Path to the base directory of the submission
    Returns:
        List[List[str]]: A list of lists containing the extracted data.
        Each inner list represents a row with the following columns:
        [sample, fastq1, fastq2, bed_file, reference]
    """
    submission_base_path = Path(submission_base_path)

    data = []
    submission = json_data["submission"]
    genomic_study_subtype = submission.get("genomicStudySubtype", "")

    donors = json_data["donors"]

    for donor in donors:
        case_id = donor.get("tanG", "")
        lab_data = donor.get("labData", [])

        for lab_datum in lab_data:
            lab_data_name = lab_datum["labDataName"]

            # derive sample_id from case_id and lab_data_name
            sample_id = f"""{case_id}_{lab_data_name.replace(" ", "_")}"""

            library_type = lab_datum["libraryType"]
            sequence_subtype = lab_datum["sequenceSubtype"]
            sequence_data = lab_datum.get("sequenceData", [])
            sequencing_layout = lab_datum["sequencingLayout"]

            files = sequence_data.get("files", [])

            # determine file structure for this lab datum
            reference = sequence_data["referenceGenome"]
            fastq_files = [f for f in files if f["fileType"] == "fastq"]

            bed_files = [f for f in files if f["fileType"] == "bed"]
            # there should be at most one bed file
            assert len(bed_files) < 2, (
                "More than one bed file specified in the submission metadata!"
            )

            if len(bed_files) > 0:
                bed_file = bed_files[0]
                bed_file_path = (
                    submission_base_path / "files" / bed_file["filePath"]
                ).absolute()
            else:
                bed_file = None
                bed_file_path = None

            if sequencing_layout == "paired-end":
                for fastq_r1, fastq_r2 in determine_fastq_pairs(fastq_files):
                    fastq_r1_file_path = (
                        submission_base_path / "files" / fastq_r1["filePath"]
                    ).absolute()
                    fastq_r2_file_path = (
                        submission_base_path / "files" / fastq_r2["filePath"]
                    ).absolute()

                    yield {
                        "sample": sample_id,
                        "labDataName": lab_data_name,
                        "libraryType": library_type,
                        "sequenceSubtype": sequence_subtype,
                        "genomicStudySubtype": genomic_study_subtype,
                        "fastq_1": str(fastq_r1_file_path),
                        "fastq_2": str(fastq_r2_file_path),
                        "bed_file": (
                            str(bed_file_path) if bed_file_path is not None else None
                        ),
                        "reference": reference,
                    }
            else:
                for fastq_file in fastq_files:
                    fastq_file_path = (
                        submission_base_path / "files" / fastq_file["filePath"]
                    ).absolute()

                    yield {
                        "sample": sample_id,
                        "labDataName": lab_data_name,
                        "libraryType": library_type,
                        "sequenceSubtype": sequence_subtype,
                        "genomicStudySubtype": genomic_study_subtype,
                        "fastq_1": str(fastq_file_path),
                        "fastq_2": None,
                        "bed_file": (
                            str(bed_file_path) if bed_file_path is not None else None
                        ),
                        "reference": reference,
                    }


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Extract data from GRZ schema JSON and create a CSV file."
    )
    parser.add_argument(
        "submission_base_path",
        help="Path to the submission base directory",
    )
    parser.add_argument(
        "--submission_metadata_json",
        help="Path to the submission metadata JSON file",
        default=None,
    )
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    submission_base_path = Path(args.submission_base_path)
    output_file = "grzqc_samplesheet.csv"
    # Set default metadata file path if not provided
    metadata_file = (
        Path(args.submission_metadata_json)
        if args.submission_metadata_json
        else submission_base_path / "metadata" / "metadata.json"
    )

    # read the metadata.json file
    try:
        with open(metadata_file, "r") as f:
            json_data = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: The file '{metadata_file}' is not a valid JSON file.")
        return
    except FileNotFoundError:
        print(f"Error: The file '{metadata_file}' was not found.")
        return

    extracted_data = list(extract_data(json_data, submission_base_path))
    samples_df = pd.DataFrame.from_records(extracted_data)

    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.max_colwidth",
        None,
    ):
        print(samples_df)

    try:
        samples_df.to_csv(output_file, index=False)
        print(f"Data has been extracted and saved to '{output_file}'")
    except PermissionError:
        print(f"Error: Permission denied when trying to write to '{output_file}'.")
    except IOError as e:
        print(
            f"Error: An I/O error occurred while writing to '{output_file}': {str(e)}"
        )


if __name__ == "__main__":
    main()
