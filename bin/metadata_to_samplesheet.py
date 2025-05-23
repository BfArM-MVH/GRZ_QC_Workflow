#!/usr/bin/env python3

import argparse
import json
from itertools import groupby
from operator import attrgetter
from pathlib import Path

from grz_pydantic_models.submission.metadata import (
    FileType,
    GrzSubmissionMetadata,
    Relation,
    SequencingLayout,
)


def sanitize(s: str) -> str:
    return s.replace(" ", "_")


def main(submission_root: Path):
    with open(submission_root / "metadata" / "metadata.json") as metadata_file:
        metadata = GrzSubmissionMetadata(**json.load(metadata_file))

    samples = []
    for donor in metadata.donors:
        targets = None

        # the tanG/VNg *CANNOT* be stored permanently
        # the donor_pseudonym for the index patient is the tanG, so redact it
        donor_pseudonym = (
            donor.donor_pseudonym if donor.relation != Relation.index_ else "index"
        )

        for lab_datum in donor.lab_data:
            sample_id = f"{donor_pseudonym}_{sanitize(lab_datum.lab_data_name)}"

            if lab_datum.sequence_data is not None:
                read_files = []

                for file in lab_datum.sequence_data.files:
                    match file.file_type:
                        case FileType.bed:
                            targets = file.file_path
                        case FileType.fastq:
                            read_files.append(file)

                if read_files:
                    if lab_datum.sequencing_layout == SequencingLayout.paired_end:
                        read_files_sorted = sorted(
                            read_files,
                            key=attrgetter("flowcell_id", "lane_id", "read_order"),
                        )
                        read_files_paired = groupby(
                            read_files_sorted, key=attrgetter("flowcell_id", "lane_id")
                        )
                        subsamples = ((r1, r2) for _meta, (r1, r2) in read_files_paired)
                    else:
                        subsamples = ((r, None) for r in read_files)

                    for reads1, reads2 in subsamples:

                        def resolve(p: Path) -> str:
                            return (
                                ""
                                if p is None
                                else str((submission_root / "files" / p).absolute())
                            )

                        samples.append(
                            {
                                "sample": sanitize(sample_id),
                                "laneId": reads1.lane_id,
                                "flowcellId": reads1.flowcell_id,
                                "labDataName": sanitize(lab_datum.lab_data_name),
                                "libraryType": sanitize(lab_datum.library_type),
                                "sequenceSubtype": sanitize(lab_datum.sequence_subtype),
                                "genomicStudySubtype": sanitize(
                                    metadata.submission.genomic_study_subtype
                                ),
                                "sequencerManufacturer": sanitize(
                                    lab_datum.sequencer_manufacturer
                                ),
                                "fastq_1": resolve(reads1.file_path),
                                "fastq_2": resolve(
                                    None if reads2 is None else reads2.file_path
                                ),
                                "bed_file": resolve(targets),
                                "reference": lab_datum.sequence_data.reference_genome,
                            }
                        )

    with open("grzqc_samplesheet.csv", "w") as output_file:
        output_file.write(
            "sample,laneId,flowcellId,labDataName,libraryType,sequenceSubtype,genomicStudySubtype,sequencerManufacturer,fastq_1,fastq_2,bed_file,reference\n"
        )
        for sample in samples:
            output_file.write(
                ",".join(
                    [
                        sample["sample"],
                        sample["laneId"],
                        sample["flowcellId"],
                        sample["labDataName"],
                        sample["libraryType"],
                        sample["sequenceSubtype"],
                        sample["genomicStudySubtype"],
                        sample["sequencerManufacturer"],
                        sample["fastq_1"],
                        sample["fastq_2"],
                        sample["bed_file"],
                        sample["reference"],
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a samplesheet from a submission"
    )
    parser.add_argument(
        "submission_root", help="Path to the submission base directory", type=Path
    )

    args = parser.parse_args()

    main(submission_root=args.submission_root)
