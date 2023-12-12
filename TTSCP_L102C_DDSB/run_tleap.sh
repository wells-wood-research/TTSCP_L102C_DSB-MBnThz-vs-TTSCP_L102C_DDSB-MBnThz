#!/bin/bash

# Iterate over all .in files in the current directory
for input_file in *.in; do
    # Check if the file exists and is a regular file
    if [ -f "$input_file" ]; then
        # Extract the base name without the extension
        base_name=$(basename -s .in "$input_file")
        
        # Run tleap with the current .in file and direct the output to a file with the same base name
        tleap -f "$input_file" > "${base_name}_output.log"
    fi
done
