# Adjust paths as needed
/Users/user/bin/rosetta_bin_linux_2020.37.61417_bundle/main/source/bin/rosetta_scripts.static.linuxgccrelease \
	-database /Users/user/bin/rosetta_bin_linux_2020.37.61417_bundle/main/database/ \
	-in:file:fullatom \
	-s $line \
	-detect_disulf true \
	-parser:protocol ./disulfide_eval.xml \
	-out:pdb true \
	-out:file:scorefile dsf_eval_${line}.sc \
	-run:preserve_header \
	-overwrite