process SUBMISSION_TO_SAMPLES_JSON {
  input:
    path submission_root

  output:
    stdout

  script:
    """
    #!/usr/bin/env python
    import itertools
    import json
    from operator import attrgetter
    from pathlib import Path

    from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata, FileType, SequencingLayout

    submission_root = Path('${submission_root}')
    with open(submission_root / 'metadata' / 'metadata.json') as metadata_file:
        metadata = GrzSubmissionMetadata(**json.load(metadata_file))

    regions = None
    samples = []

    for donor in metadata.donors:
        for i, lab_datum in enumerate(donor.lab_data):
            if lab_datum.sequence_data is not None:
                read_files = []
                read_file_types = []

                for file in lab_datum.sequence_data.files:
                    match file.file_type:
                        case FileType.bed:
                            regions = file.file_path
                        case FileType.bam | FileType.fastq:
                            read_files.append(file)
                            read_file_types.append(file.file_type)

                if read_files:
                    if lab_datum.sequencing_layout == SequencingLayout.paired_end:
                        is_single_end = False
                        read_files_sorted = sorted(read_files, key = attrgetter('flowcell_id', 'lane_id', 'read_order'))
                        read_files_paired = itertools.groupby(read_files_sorted, key = attrgetter('flowcell_id', 'lane_id'))
                        read_paths = list(itertools.chain.from_iterable(((r1.file_path, r2.file_path) for k, (r1, r2) in read_files_paired)))
                    else:
                        is_single_end = True
                        read_paths = [file.file_path for file in read_files]

                    sample = {
                        'meta': {
                            'id': f"{donor.relation}{i}",
                            'donor_relation': donor.relation,
                            'lab_data_name': lab_datum.lab_data_name,
                            'library_type': lab_datum.library_type,
                            'sequence_subtype': lab_datum.sequence_subtype,
                            'reference_genome': lab_datum.sequence_data.reference_genome,
                            'single_end': is_single_end,
                            'genomic_study_subtype': metadata.submission.genomic_study_subtype
                        },
                        'reads': read_paths,
                        'file_types': read_file_types,
                        'regions': regions
                    }
                    samples.append(sample)

    print(json.dumps(samples))
    """
}

workflow PARSE_SUBMISSION {
  take:
    submission_root

  emit:
    samples = SUBMISSION_TO_SAMPLES_JSON(submission_root)
      .splitJson()
      .map{ sample ->
        [
          "meta": sample.meta,
          "reads": sample.reads.collect{ read_path -> submission_root / 'files' / read_path },
          "read_file_types": sample.file_types,
          "regions": sample.regions ? submission_root / 'files' / sample.regions : null
        ]
      }
}
