# Adjust source paths as needed
# Adjust prefix as needed
/Users/user/rosetta_bin_linux_2020.37.61417_bundle/main/source/bin/rosetta_scripts.static.linuxgccrelease \
	-database /Users/user/bin/rosetta_bin_linux_2020.37.61417_bundle/main/database/ \
	-in:file:fullatom \
	-s ../peptide_template.pdb \
	-detect_disulf true \
	-parser:protocol ../Topobuild_dsf_staple.xml \
	-out:pdb true \
	-out:prefix L2E3L2E4L3H8L3_${array}_ \
	-out:file:scorefile L2E3L2E4L3H8L3_${array}_prelim_topo.sc \
	-run:preserve_header \
	-overwrite