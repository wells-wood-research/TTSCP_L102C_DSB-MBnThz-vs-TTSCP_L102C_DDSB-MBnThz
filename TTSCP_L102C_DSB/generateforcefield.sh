#!/bin/bash

# Check if parmchk2 is available in the PATH
if ! command -v parmchk2 &> /dev/null; then
    echo "Error: parmchk2 is not found in the PATH. Please make sure AmberTools is installed and sourced."
    exit 1
fi

# Loop through all .mol2 files in the current directory
for mol2_file in *.mol2; do
    # Extract the base filename (without the extension)
    base_filename="${mol2_file%.mol2}"

    # Check if a corresponding .frcmod file already exists
    frcmod_file="${base_filename}.frcmod"
    if [ -e "$frcmod_file" ]; then
        echo "Skipping $mol2_file: $frcmod_file already exists."
    else
        # Run parmchk2 to generate the .frcmod file
        parmchk2 -i "$mol2_file" -f mol2 -o "$frcmod_file" -a Y
        if [ $? -eq 0 ]; then
            echo "Generated $frcmod_file from $mol2_file."
        else
            echo "Error generating $frcmod_file from $mol2_file."
        fi
    fi
done

