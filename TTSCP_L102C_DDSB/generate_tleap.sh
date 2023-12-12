#!/bin/bash

# Assuming mol2 and frcmod files have consistent names like MBNcage.mol2 and MBNcage.frcmod
for base_name in $(ls -1 *.mol2 | sed 's/.mol2//'); do
    input_file="${base_name}.in"
    cat <<EOF > "$input_file"
source leaprc.protein.ff14SB
source leaprc.gaff
${base_name} = loadmol2 ${base_name}.mol2
check ${base_name}
loadamberparams ${base_name}.frcmod
check ${base_name}
saveoff ${base_name} ${base_name}.lib
saveamberparm ${base_name} ${base_name}.prmtop ${base_name}.rst7
quit
EOF
done

