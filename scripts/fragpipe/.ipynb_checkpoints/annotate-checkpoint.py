#!/usr/bin/env python3

'''
This script is used to replace Accession Numbers with Annotations in the column `Annotations`. The script can only be used on the formatted proteome files that have the extension *_formatted.txt.

USAGE: python annotate.py </path/annotations_file.txt> </path/input_foldername/>
'''

import glob
import os
import sys

# Reading in the SPO numbers and associated annotations and saving them in the dictionary SPO_dict
input_file = sys.argv[1]
SPO_dict = {}

with open(input_file, 'r') as file:
    for line in file:
        spo, annotation = line.strip().split('\t', 1)  # Split the line at the first tab character
        SPO_dict[spo] = annotation

# Function to annotate the files (Replacing values in column 2) if annotations are not provided
def annotate(folder):
    for filepath in glob.glob(os.path.join(folder, '*_formatted.txt')):
        filename = os.path.basename(filepath)  # Get the filename
        print(f"Processing file: {filename}")
    
        output_file_path = filepath.replace('_formatted.txt', '_annotated.txt')
    
        with open(filepath, 'r') as input_file, open(output_file_path, 'w') as output_file:
            # Write the header line to the output file
            header = input_file.readline().strip() + '\n'
            output_file.write(header)
        
            # Process each line in the input file
            for line in input_file:
                # Split the line into columns
                columns = line.strip().split('\t')
                spo_id = columns[0]  # Extract the SPO ID from column 1
            
                # Look up the annotation in the dictionary
                annotation = SPO_dict.get(spo_id, 'Unknown')
                columns[1] = annotation
            
                # Join the columns back into a line and write to the output file
                output_line = '\t'.join(columns) + '\n'
                output_file.write(output_line)
    
        print(f"Annotated lines written to: {output_file_path}")


if __name__ == "__main__": 
    
    if len(sys.argv) !=3: # If there is more than 2 argument to call this script, it will provide guidance on how to use it.
        print("Usage: python annotate.py </path/annotations_file.txt> </path/input_foldername/>")
        sys.exit(1)
    
    folder = sys.argv[2]
    annotate(folder)