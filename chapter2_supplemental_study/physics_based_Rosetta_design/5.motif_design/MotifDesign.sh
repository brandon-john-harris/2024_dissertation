#!/bin/bash
# adjust paths as needed
ROSETTA_PATH=/Users/user/rosetta_bin_linux_2020.37.61417_bundle/main
${ROSETTA_PATH}/source/bin/rosetta_scripts.static.linuxgccrelease \
	-database ${ROSETTA_PATH}/database/ \
	-in:file:fullatom \
	-nstruct 1 \
	-l design_inputs.txt \
	-detect_disulf true \
	-parser:protocol ./seed_files/MotifDesign.xml \
	-out:pdb true \
	-no_optH false \
	-ignore_unrecognized_res \
	-out:prefix ${array_prefix}_MotifDesign_ \
	-run:preserve_header \
	-out:level 300 \
	-overwrite
