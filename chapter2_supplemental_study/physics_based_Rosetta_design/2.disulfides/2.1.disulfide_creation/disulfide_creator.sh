# Adjust paths as needed
/Users/user/bin/rosetta_bin_linux_2020.37.61417_bundle/main/source/bin/rosetta_scripts.static.linuxgccrelease \
	-database /Users/user/bin/bin/rosetta_bin_linux_2020.37.61417_bundle/main/database/ \
	-in:file:fullatom \
	-s $line \
	-detect_disulf true \
	-parser:protocol ./disulfide_create.xml \
	-out:pdb true \
	-run:preserve_header \
	-overwrite
