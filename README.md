# Scripts for building the 3D ex vivo MRI atlas 

## Scripts

1. "copy_inputs_flywheel" downdloads the scans from flywheel and does preprocessing. Perform n4 correction first. Then create manual raw-to-axisalign affine transformation. This script contains the function to apply the affine transform to n4 corrected scans. 

2. "organize_inputs" organizes the data required for atals construction in a preproc directory

3. "build_mst_srlm_atlas2022" is the main script for constructing the atlas. 

4. "reslice_histo_segs" maps the histology annotations from subject space to the template and computes the consensus segmentation

5. tau_density_to_atlas maps the NFT burden maps from subject space to the templae and computes the average maps and ROI measures

6. template_skel_analysis contains the thickness computation and meshglm analyses for surface based analyses

 	





