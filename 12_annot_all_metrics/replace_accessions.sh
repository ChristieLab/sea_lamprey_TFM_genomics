#!/bin/bash

# This script replaces FASTA headers in an input file with
# corresponding parent RNA accessions from a lookup file.

# Check if the correct number of arguments is provided.
if [ "$#" -ne 3 ]; then
    echo "Usage: ./replace_fasta_headers.sh <input_sequence_file> <parent_accessions_file> <output_file>"
    exit 1
fi

input_sequence_file="$1"
parent_accessions_file="$2"
output_file="$3"

# Check if input files exist.
if [ ! -f "$input_sequence_file" ]; then
    echo "Error: Input sequence file '$input_sequence_file' not found."
    exit 1
fi

if [ ! -f "$parent_accessions_file" ]; then
    echo "Error: Parent accessions file '$parent_accessions_file' not found."
    exit 1
fi

# Use awk to process both files in a single pass.
awk '
    # Process the parent file first to build a lookup table.
    # FNR is the record number in the current file, NR is the overall record number.
    FNR==NR {
        # The parent accession is after "rna-".
        # We capture the part that will be used for matching.
        match($0, /rna-(XM_[0-9]+\.[0-9]+)/, arr);
        if (arr[1] != "") {
            # Store the full parent name with the core accession as the key.
            parents[arr[1]] = $0;
        }
        next;
    }

    # Process the sequence file.
    /^>/ {
        # If the line is a FASTA header, extract the core accession.
        match($0, /(XM_[0-9]+\.[0-9]+)/, arr);
        accession = arr[1];

        # Look up the core accession in the parent array.
        if (accession in parents) {
            # If a match is found, print the new header.
            print ">" parents[accession];
        } else {
            # If no match, print the original header.
            print $0;
        }
    }

    !/^>/ {
        # If the line is not a header (it is sequence data), print it as is.
        print $0;
    }
' "$parent_accessions_file" "$input_sequence_file" > "$output_file"

echo "Successfully replaced headers. New file saved as '$output_file'."