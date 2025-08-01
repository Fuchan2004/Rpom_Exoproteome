#!/usr/bin/env python3

'''
This script is used to format fragpipe output files such that they are suitable for running in with the scripts 'statistics.py', 'xy-plots.py', and 'volcano_plot.py'. 
You obtain a file called 'report.pg_matrix.tsv' file from fragpipe containing the normalized DIA data. The column names are as follows: 
Protein.Group / Protein.Ids / Protein.Names / Genes / Frist.Protein.Description / Sample 1 / Sample 2 / etc. However, for the below script to work, make sure that: 
    1) there is only one sample (in triplicates) present in your file. Thus the file has 8 columns with the last three being Sample_1_1 / Sample_1_2 / Sample_1_3 .
    2) the file name has following format: 'DATE_STRAIN_MEDIUM_GROWTHPHASE.txt'
    3) the annotations are listed under Protein.Ids. If there are no annotions but accession numbers instead, use the script 'annotate.py' to replace accession numbers with the annotations. 

USAGE: python format_fragpipe_output.py <foldername>
'''


import glob
import os
import sys


def format_proteomefiles(folder=''):
    
    # Define the column indices to keep: [Accession, Annotation, triplicates]
    selected_columns = [0, 1, 5, 6, 7]
    
    for filepath in glob.glob(os.path.join(folder, '*.txt')):
        filename = os.path.basename(filepath)
        print(f"Processing file: {filename}")

        # Read the original file
        with open(filepath, 'r') as file:
            lines = file.readlines()

        # Rename columns to match downstream script expectations
        header = lines[0].strip().split('\t')

        # Ensure there are at least 8 columns (Accession + triplicates)
        if len(header) == 8:
            header[0] = 'Accession Number'
            header[1] = 'Annotation'
            header[5] = '1'
            header[6] = '2'
            header[7] = '3'
        else:
            print(f"Skipping {filename}: Incorrect format.")
            continue

        # Filter selected columns for the header
        filtered_header = [header[i] for i in selected_columns]
        formatted_lines = ['\t'.join(filtered_header) + '\n']

        # Read and filter the data rows
        for line in lines[1:]:
            columns = line.strip().split('\t')
            
            # Ensure row has at least enough columns to access the selected ones
            filtered_columns = []
            for idx in selected_columns:
                if idx < len(columns):
                    filtered_columns.append(columns[idx])
                else:
                    # Add empty string if column is missing
                    filtered_columns.append('')

            formatted_lines.append('\t'.join(filtered_columns) + '\n')

        # Sort by Accession Number (first column)
        sorted_lines = [formatted_lines[0]] + sorted(formatted_lines[1:], key=lambda x: x.split('\t')[0])

        # Write to a new file with "_formatted" suffix
        output_file_path = filepath.replace('.txt', '_formatted.txt')
        
        with open(output_file_path, 'w') as file:
            file.writelines(sorted_lines)

        print(f"Formatted file written to: {output_file_path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python format_fragpipe_output.py <input_foldername>")
        sys.exit(1)

    folder = sys.argv[1]
    format_proteomefiles(folder)