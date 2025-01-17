#!/usr/bin/env python

import json
import csv
import argparse
import re
import os
from typing import List, Dict

def extract_data(json_data: Dict, submission_base_path) -> List[List[str]]:
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
    data = []
    donors = json_data['donors']

    for donor in donors:
        
        case_id = donor.get('caseId', '')
        print(case_id)
        lab_data_list = donor.get('labData', [])

        for lab_data in lab_data_list:
            #print("----------Lab data---------------")
            #print(lab_data)
            lab_data_name = lab_data.get('labDataName', '').replace(' ', '_')
            libraryType = lab_data.get('libraryType', [])
            sequence_data = lab_data.get('sequenceData', [])
            
            #print(sequence_data)
            #for sequence_data in sequence_data_list:
            #print("----------Sequence data---------------")
            #print(sequence_data)
            files = sequence_data.get('files', [])
            sample = f"{case_id}_{lab_data_name}_{libraryType}"
            bed_file = ""
            reference= sequence_data.get('referenceGenome', [])
            #print(reference)
            fastq_files = {"R1": [], "R2": []}
            for file in files:
                file_path = submission_base_path+"/files/"+file.get('filePath', '')
                file_type = file.get('fileType', '')

                if file_type == 'fastq':
                    if re.search(r'(R1|read1)', file_path, re.IGNORECASE):
                        fastq_files["R1"].append(file_path)
                    elif re.search(r'(R2|read2)', file_path, re.IGNORECASE):
                        fastq_files["R2"].append(file_path)
                elif file_type == 'bed':
                    bed_file = file_path

            # Match R1 and R2 files
            for r1_file in fastq_files["R1"]:
                r1_base = re.sub(r'(R1|read1)', '', r1_file, flags=re.IGNORECASE)
                for r2_file in fastq_files["R2"]:
                    r2_base = re.sub(r'(R2|read2)', '', r2_file, flags=re.IGNORECASE)
                    if r1_base == r2_base:
                        data.append([sample, r1_file, r2_file, bed_file, reference])
                        break
    return data

def main():
    parser = argparse.ArgumentParser(description='Extract data from GRZ schema JSON and create a CSV file.')
    parser.add_argument('submission_base_path', help='Complete path to the root directory of the submission')
    parser.add_argument('output_path', help='Complete path to the root directory of the submission')
    args = parser.parse_args()

    # Check if the submission base path exists
    if not os.path.isdir(args.submission_base_path):
        print(f"Error: The submission base path '{args.submission_base_path}' does not exist.")
        return

    metadata_file = args.submission_base_path+"metadata/metadata.json"
    output_file = args.output_path+"grzqc_samplesheet.csv"
    #print(metadata_file)
    
    try:
        with open(metadata_file, 'r') as f:
            json_data = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: The file '{metadata_file}' is not a valid JSON file.")
        return
    except FileNotFoundError:
        print(f"Error: The file '{metadata_file}' was not found.")
        return
    
    

    #print(json_data['submission']['submissionDate'])
    extracted_data = extract_data(json_data, args.submission_base_path)
    #print(extracted_data)
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['sample', 'fastq_1', 'fastq_2', 'bed_file', 'reference'])
            writer.writerows(extracted_data)
        print(f"Data has been extracted and saved to '{output_file}'")
    except PermissionError:
        print(f"Error: Permission denied when trying to write to '{output_file}'.")
    except IOError as e:
        print(f"Error: An I/O error occurred while writing to '{output_file}': {str(e)}")

if __name__ == "__main__":
    main()
