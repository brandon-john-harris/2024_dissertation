#!/bin/bash
# adjust paths as needed
# Updated 2021-05-05 Script created

if [ $# -lt 2 ]; then
    echo "USAGE: ./positive_score_filter.sh <PDB_flist_file> <array #>"
    exit
fi
pdb_file=${1}
array_number=${2}
/Users/user/rosetta_bin_linux_2020.37.61417_bundle/main/source/bin/rosetta_scripts.default.linuxgccrelease \
	-database /Users/user/rosetta_bin_linux_2020.37.61417_bundle/main/database/ \
	-in:file:fullatom \
	-l $pdb_file \
	-detect_disulf true \
	-parser:protocol ./score_filter.xml \
	-out:pdb true \
	-out:file:scorefile positive_score_filter_${array_number}.sc \
	-run:preserve_header \
	-overwrite