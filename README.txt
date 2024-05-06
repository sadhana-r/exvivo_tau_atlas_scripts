
*copy_inputs_flywheel.sh	Created in Jan 2022. Downloads raw scans from flywheel. Perform n4 correction first. Then create manual raw-to-axisalign affine transformation. And contains function to apply the affine transform to n4 corrected scans. 

*build_mst: 	atlas building script used to build initial atlas (24 subjects: subj24_08072019.txt) for first paper submission (May 2020). This version does not include the SRLM label

*build_mst_srlm:	Atlas building script used for Acta Neuro Comm 2021 paper (29 subjects - subj_29.txt). Includes the srlm label in the registration pipeline.
	* add_subj_to_template.sh	Built on this script by adding new subjects to the template for AAIC 2022 abstract. List of subjects in this dataset:subj_aaic2022.txt




