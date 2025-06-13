#!/usr/bin/env python3

import argparse
import json
from itertools import groupby
from operator import attrgetter
from pathlib import Path

import pandas as pd
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
    for donor_index, donor in enumerate(metadata.donors):
        targets = None

        # the tanG/VNg *CANNOT* be stored permanently
        # the donor_pseudonym for the index patient is the tanG, so redact it
        donor_pseudonym = (
            donor.donor_pseudonym if donor.relation != Relation.index_ else "index"
        )

        for lab_datum_index, lab_datum in enumerate(donor.lab_data):
            sample_id = f"{donor.relation}{donor_index}_{lab_datum.sequence_subtype}{lab_datum_index}"

            if lab_datum.sequence_data is not None:
                read_files = []

                for file in lab_datum.sequence_data.files:
                    match file.file_type:
                        case FileType.bed:
                            targets = file.file_path
                        case FileType.fastq | FileType.bam:
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
                        runs = ((r1, r2) for _meta, (r1, r2) in read_files_paired)
                    else:
                        runs = ((r, None) for r in read_files)

                    for run_index, (reads1, reads2) in enumerate(runs):
                        run_id = ""

                        if reads1.lane_id:
                            run_id = sanitize(reads1.lane_id)

                        if reads1.flowcell_id:
                            flowcell_id_san = sanitize(reads1.flowcell_id)
                            run_id = (
                                run_id + f"_{flowcell_id_san}"
                                if run_id
                                else flowcell_id_san
                            )

                        def resolve(p: Path) -> str:
                            return (
                                ""
                                if p is None
                                else str((submission_root / "files" / p).absolute())
                            )

                        samples.append(
                            {
                                "sample": sample_id,
                                "runId": run_id if run_id else f"run{run_index}",
                                "donorPseudonym": donor_pseudonym,
                                "laneId": reads1.lane_id,
                                "flowcellId": reads1.flowcell_id,
                                "labDataName": lab_datum.lab_data_name,
                                "libraryType": lab_datum.library_type,
                                "sequenceSubtype": lab_datum.sequence_subtype,
                                "genomicStudySubtype": metadata.submission.genomic_study_subtype,
                                "sequencerManufacturer": lab_datum.sequencer_manufacturer,
                                "reads_long": resolve(reads1.file_path)
                                if lab_datum.library_type.endswith("_lr")
                                else None,
                                "reads1": resolve(reads1.file_path)
                                if not lab_datum.library_type.endswith("_lr")
                                else None,
                                "reads2": resolve(
                                    None if reads2 is None else reads2.file_path
                                )
                                if not lab_datum.library_type.endswith("_lr")
                                else None,
                                "bed_file": resolve(targets),
                                "reference": lab_datum.sequence_data.reference_genome,
                                "meanDepthOfCoverage": lab_datum.sequence_data.mean_depth_of_coverage,
                                "targetedRegionsAboveMinCoverage": lab_datum.sequence_data.targeted_regions_above_min_coverage,
                                "percentBasesAboveQualityThreshold": lab_datum.sequence_data.percent_bases_above_quality_threshold.percent,
                            }
                        )

    samples_df = pd.DataFrame(samples)
    samples_df.to_csv("grzqc_samplesheet.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a samplesheet from a submission"
    )
    parser.add_argument(
        "submission_root", help="Path to the submission base directory", type=Path
    )

    args = parser.parse_args()

    main(submission_root=args.submission_root)
