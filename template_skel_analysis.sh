#/bin/bash
set -x -e

module load R/3.2.5

ROOT=/project/hippogang_3/sravikumar/atlasPHG2019
ADIR=$ROOT/atlasPHG_V3_2022
WORK=$ADIR/shape_animate
THKDIR=$WORK/thickness

# Read the common include file
source $ROOT/scripts/common.sh

SF_NAMES=("NULL" "CA1" "CA2" "DG" "CA3" "HATA" "SUB" "CYST" "MISC" "SRLM")


#TEMPLATEDIR=$ADIR/mst_multilabel/template_ncc/iter5
#TEMPLATEDIR=$ADIR/mst_multilabel_v2/template_ncc/iter5
TEMPLATEDIR=$ADIR/mst_multilabel/template_ncc/iter5

#TEMPLATE_SEG=$TEMPLATEDIR/template_phg_snake_edited.nii.gz
TEMPLATE_SEG=$TEMPLATEDIR/template_phg_thresh0.5.nii.gz
TEMPLATE_SRLM=$TEMPLATEDIR/template_srlm_thresh0.5.nii.gz

SUBJ_LIST=$ROOT/scripts/subj_atlas_2022_v2.txt

mkdir -p $WORK

# Common steps for the thickness analysis
function thickness_templates()
{
  # Extract the surface of the dark band for thickness analyses
  mkdir -p $THKDIR
 

 # Extract the surface of the srlm for thickness analysis
  vtklevelset $TEMPLATE_SRLM $THKDIR/template_srlm.vtk 0.5 
  
 #Extract the surface of the gray matter - srlm for thickness analyses
  c3d $TEMPLATE_SRLM -thresh 0.5 inf 1 0 \
  $TEMPLATE_SEG -add -retain-labels 1 $TMPDIR/graystrip.nii.gz

  vtklevelset $TMPDIR/graystrip.nii.gz $THKDIR/template_gm.vtk 0.5


  # Apply Taubin smoothing to the meshes and create the skeleton mesh
   for what in srlm; do #gm

   	/project/hippogang_2/pauly/bin/mesh_smooth_curv -mu -0.51 -lambda 0.5 -iter 1000 -r 2.0 \
  	$THKDIR/template_${what}.vtk \
  	$THKDIR/template_${what}_tsm.vtk
   
    # Strip point data that messes Greedy up
  	cat $THKDIR/template_${what}_tsm.vtk \
  	| awk '$1 == "POINT_DATA" {stop=1} { if (stop!=1) print $0  }' \
  	> $THKDIR/template_${what}_tsm_nopd.vtk

    #Reduce pruning to 1.2  - thickness skeleton is used for sampling from the thickness images
    #Try set pruning to 2 - April 5 2023
    $CMREP_LMSHOOT/cmrep_vskel \
    -Q /project/hippogang_3/sravikumar/atlasPHG2019/pkgs/libc/qvoronoi \
    -e 5 -p 2.0 -c 1 \
    $THKDIR/template_${what}_tsm_nopd.vtk $THKDIR/template_${what}_tsm_skel.vtk


  if [[ $what == "gm" ]]; then
    #Probability of each voxel's label is based on its distance to that region. -dup -times -sqrt computes the absolute value of the distance
    c3d $ADIR/histology_annotation_coarse/template_hf_label.nii.gz -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz

    $ROOT/pkgs/mesh_image_sample -B -i 0 $THKDIR/template_gm_tsm_skel.vtk $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz $THKDIR/template_gm_tsm_skel_labeled.vtk HF_label

    #Full histology annotation
    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -replace 5 0 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $THKDIR/template_gm_tsm_skel_labeled.vtk $TMPDIR/seg_expanded_coarse.nii.gz $THKDIR/template_gm_tsm_skel_labeled.vtk histo_label_coarse

    #Fine histo seg
    c3d -verbose $ADIR/histology_annotation_fine/template_histo_seg_final.nii.gz -replace 5 0 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_fine.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $THKDIR/template_gm_tsm_skel_labeled.vtk $TMPDIR/seg_expanded_fine.nii.gz $THKDIR/template_gm_tsm_skel_labeled.vtk histo_label_fine

  elif [[ $what == "srlm" ]]; then

    #Full histology annotation
    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -retain-labels 5 13 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $THKDIR/template_srlm_tsm_skel.vtk $TMPDIR/seg_expanded_coarse.nii.gz $THKDIR/template_srlm_tsm_skel_labeled.vtk histo_label_coarse

  fi

  done
}

# Subject-specific ROI thickness analysis - extract mean and median ROI thickness
# Do this using image space thickness maps and segmentations - Instead of dumping the skeleton mesh
function thickness_perlabel()
{
  id=${1?}
  what=${2?}

  #VTK library needed for mesh_fill_missing_data
  export LD_LIBRARY_PATH=/project/hippogang_3/sravikumar/packages/vtk-release-6.3/bin/lib

  # Subject's final directory
  SDIR=$ADIR/mst_multilabel/final/${id}

  # Thickness per subject directory
  WDIR=$THKDIR/subj_skel/${id}

  #Output directory
  ODIR=$THKDIR/thickness_ROI
  mkdir -p $ODIR

  #Output segmentation in subject space - template warped to subject space. SRLM consists of label 5 and 13 (PP)
  if [[ $what == "gm" ]]; then
    SUBJSEG="$WORK/template_tosubj/${id}/${id}_unwarp_histoseg_coarse.nii.gz -replace 5 0 13 0"
  elif [[ $what == "srlm" ]]; then
    SUBJSEG="$WORK/template_tosubj/${id}/${id}_unwarp_histoseg_coarse.nii.gz -retain-labels 5 13"
  fi

  #Thickness image
  local THKIMG_MD=$WDIR/${id}_${what}_thickness_md.mha

  #Instead of using subject space segmentation with missing data label (from tearing and bad registration),
  # Use the thickness image thresholded since there may be more missing data points from tetfill.

  # Giving image region errors!! Image origins are 
  # # very slightly off!!!!!!! why!!!! Need to use copy-transform to make sure origins are exactly the same (hack)
  c3d $SUBJSEG -interpolation NearestNeighbor -resample 300% -as HSEG \
  $THKIMG_MD -thresh 0.001 inf 1 0 -as THKMASK \
  -push HSEG -reslice-identity -popas HRSL -push THKMASK -push HRSL -copy-transform \
  -push THKMASK -times -o $TMPDIR/${id}_${what}_histoseg_masked.nii.gz

  c3d $THKIMG_MD $TMPDIR/${id}_${what}_histoseg_masked.nii.gz -lstat > $ODIR/${id}_${what}_lstat_thickness_coarse.txt
}


# Subject-specific ROI thickness analysis - extract mean and median ROI thickness
# Do this using image space thickness maps and segmentations - Instead of dumping the skeleton mesh
function thickness_perlabel_EC()
{
  id=${1?}

  #VTK library needed for mesh_fill_missing_data
  export LD_LIBRARY_PATH=/project/hippogang_3/sravikumar/packages/vtk-release-6.3/bin/lib

  # Subject's final directory
  SDIR=$ADIR/mst_multilabel/final/${id}

  # Thickness per subject directory
  WDIR=$THKDIR/subj_skel/${id}

  #Output directory
  ODIR=$THKDIR/thickness_ROI
  mkdir -p $ODIR

  #Output segmentation in subject space - template warped to subject space. Only interested in ERC subfields
  SUBJSEG="$WORK/template_tosubj/${id}/${id}_unwarp_histoseg_fine.nii.gz -retain-labels 27 30 31 32 32 33 34 35 36"

  #Thickness image
  local THKIMG_MD=$WDIR/${id}_gm_thickness_md.mha

  #Instead of using subject space segmentation with missing data label (from tearing and bad registration),
  # Use the thickness image thresholded since there may be more missing data points from tetfill.

  # Giving image region errors!! Image origins are 
  # # very slightly off!!!!!!! why!!!! Need to use copy-transform to make sure origins are exactly the same (hack)
  c3d $SUBJSEG -interpolation NearestNeighbor -resample 300% -as HSEG \
  $THKIMG_MD -thresh 0.001 inf 1 0 -as THKMASK \
  -push HSEG -reslice-identity -popas HRSL -push THKMASK -push HRSL -copy-transform \
  -push THKMASK -times -o $TMPDIR/${id}_gm_histoseg_masked.nii.gz

  c3d $THKIMG_MD $TMPDIR/${id}_gm_histoseg_masked.nii.gz -lstat > $ODIR/${id}_gm_lstat_thickness_ECsubfields.txt
}


# Subject-specific ROI thickness analysis - extract mean and median ROI thickness
# Uses skeleton mesh generated from thicknes_persubj step
function thickness_perlabel_badregmask()
{

  id=${1?}
  what=${2?}

  #VTK library needed for mesh_fill_missing_data
  export LD_LIBRARY_PATH=/project/hippogang_3/sravikumar/packages/vtk-release-6.3/bin/lib

  # Subject's final directory
  SDIR=$ADIR/mst_multilabel/final/${id}

   # Output directory
  WDIR=$THKDIR/subj/${id}
  mkdir -p $WDIR/reg_masked

  #Output segmentation in subject space
  SUBJSEG=$WORK/template_tosubj/${id}/${id}_unwarp_histoseg_coarse.nii.gz

  # The output skeleton mesh
  local LTHK=$WDIR/reg_masked/${id}_${what}_histo_labthick_md.txt
  local THKMESH=$WDIR/${id}_${what}_thickness.vtk

 local LABLIST
 if [[ $what == "gm" ]]; then
        LABLIST="6 11 13 14 20 21 22 23 26 27 37 42 45 46 51 53 59"
	local THKMESH_MD=$WDIR/${id}_${what}_thickness_md_badreg.vtk
 elif [[ $what == "srlm" ]]; then
        LABLIST="5"
	local THKMESH_MD=$WDIR/${id}_${what}_thickness_md.vtk
 fi
   
  # Write all the labels into a file
   dumpmeshattr $THKMESH_MD Thickness_nan > $TMPDIR/${id}_${what}_dumpthick.txt
   dumpmeshattr $THKMESH_MD HISTO_LABEL > $TMPDIR/${id}_${what}_dumplabel.txt
   dumpmeshattr $THKMESH AreaElement > $TMPDIR/${id}_${what}_dumparea.txt
   paste \
	$TMPDIR/${id}_${what}_dumpthick.txt \
        $TMPDIR/${id}_${what}_dumplabel.txt \
        $TMPDIR/${id}_${what}_dumparea.txt \
        > $LTHK

   Rscript $ROOT/scripts/calculate_label_thk.R $LTHK ${id} > $WDIR/reg_masked/${id}_${what}_thickness_by_histolabel_badreg.txt

}

# Subject-specific steps for the thickness analysis
function thickness_persubj_skel()
{
  export LD_LIBRARY_PATH=/project/hippogang_3/sravikumar/packages/vtk-release-6.3/bin/lib

  id=${1?}
  what=${2?}

  # Segmentation used for this
  SEG=$TEMPLATE_SEG

  #Template skeleton
  TSKELETON=$THKDIR/template_${what}_tsm_skel.vtk

  # Subject's final directory
  SDIR=$ADIR/mst_multilabel/final/${id}

  # Output directory
  WDIR=$THKDIR/subj_skel/${id}

  # Make the directory
  mkdir -p $WDIR

  # Pad the imagebefore generating boundary - close to the edge boundary was getting cut off. Padding now happens in another function
  # The subject's segmentations - use the segmentation with the green label (to exclude from thickness)
  SUBJ_GM=$THKDIR/subj/${id}/${id}_unwarp_missinglabel_withbadreg.nii.gz

  #Use padded version for all cases
  greedy -d 3 \
    -rf $THKDIR/subj/${id}/${id}_unwarp_img_padded.nii.gz \
    -ri LABEL 0.24vox \
    -rm $ROOT/preproc/$id/${id}_axisalign_srlmseg_sr.nii.gz $THKDIR/subj/${id}/${id}_unwarp_srlmseg_sr_padded.nii.gz \
    -r $ROOT/preproc/$id/unwarp/${id}_unwarp.mat $ROOT/preproc/$id/${id}_transform_to_axisalign.mat,-1

  SUBJ_SRLM=$THKDIR/subj/${id}/${id}_unwarp_srlmseg_sr_padded.nii.gz

  # 1 - GM, 2 - missing regions, 3 - SRLM
  #Dilate and erode the srlm segmentation to fill holes in the segmentation.
  #Not an issue when warping the template to subject space  
  c3d $SUBJ_SRLM -dilate 1 5x5x5vox -erode 1 5x5x5vox -as SRLM \
  -o $TMPDIR/${id}_${what}_srlmseg_missinglabel.nii.gz \
  -thresh 1 inf 1 0 -o $WDIR/srlm_seg.nii.gz 

  # Set the name of the output segmentation file used in the next step. SRLM - 3. MD (SRLM)-4, GM - 1, MD (GM) - 2 #
  SUBJSEG=$WDIR/${id}_unwarp_phg_wsrlm.nii.gz

  #srlm_seg and gm_seg are binary images (just 1/0 - no missing data)
  #April 2023 - add -comp to make sure gm seg has one component
  c3d $WDIR/srlm_seg.nii.gz -scale 2 $SUBJ_GM -add -o $SUBJSEG -thresh 1 2 1 0 \
  -as GM -comp -thresh 1 1 1 0 -push GM -times -o $WDIR/gm_seg.nii.gz

  # The output skeleton mesh
  local RESLICED_SKEL=$WDIR/${id}_${what}_warped_skeleton.vtk
  local SKEL=$WDIR/${id}_${what}_skeleton.vtk
  local THKIMG=$WDIR/${id}_${what}_thickness.nii.gz
  local DEPTHIMG=$WDIR/${id}_${what}_depth.nii.gz
  local SUBJ_TETRA=$WDIR/${id}_${what}_tetra.vtk
  local LTHK=/tmp/labthick.txt
  local THKIMG_MASKED=$WDIR/${id}_${what}_thickness_md.mha
  local THK_SKEL=$WDIR/${id}_${what}_warped_thickness_skeleton.vtk

  # Map the dark band and gray matter skeletons from template space into subject space- SHOULD USE TSM
  $GREEDY_HOME/greedy -d 3 \
  -rf $SEG \
  -rs $THKDIR/template_${what}_tsm_skel.vtk $RESLICED_SKEL \
  -r $(cat $SDIR/chain_unwarped_to_template.txt)

  #Extract the surface from the subject segmentation
  vtklevelset $WDIR/${what}_seg.nii.gz $TMPDIR/${what}_subj.vtk 0.5

  mesh_smooth_curv -mu -0.51 -lambda 0.5 -iter 1000 \
  -r 2.0 $TMPDIR/${what}_subj.vtk $WDIR/${what}_subj_smooth.vtk


 # Extract the skeleton and compute thickness map
  if [[ ! -f $SKEL ]]; then 
    $CMREP_LMSHOOT/cmrep_vskel -Q /project/hippogang_3/sravikumar/atlasPHG2019/pkgs/libc/qvoronoi \
    -e 5 -p 1.2 -d $SUBJ_TETRA \
    -I $WDIR/${what}_seg.nii.gz $THKIMG $DEPTHIMG \
    $WDIR/${what}_subj_smooth.vtk $SKEL 
 fi

 if [[ ${what} == "srlm" ]]; then
  SUBJSEG=$TMPDIR/${id}_${what}_srlmseg_missinglabel.nii.gz
 fi

 ## As of April 2023, sample thickness on the resliced template skeleton directly from the tetrahedral instead of 
 ## mapping the thickness measures from the tetmesh to image space first
 #Instead of -d 0.4
 mesh_tetra_sample -B -d 0.6 -D "DistanceToTet" $RESLICED_SKEL $SUBJ_TETRA $THK_SKEL "VoronoiRadius"

# # Map the missing label to the skeleton to make where the thickness should be ignored
c3d -verbose $SUBJSEG -retain-labels 1 2 -trim 5vox -replace 0 1000 \
-split \
-foreach -sdt -dup -times -sqrt -reciprocal -endfor \
-scale 0 -merge \
-replace 1000 0 -o $WDIR/seg_expanded.nii.gz

# Sample the mesh by the segmentation
$ROOT/pkgs/mesh_image_sample -B -i 0 $THK_SKEL $WDIR/seg_expanded.nii.gz $THK_SKEL MissingData

python AddMissingSkelThickness.py -i $THK_SKEL -o $THK_SKEL -d "VoronoiRadius" -m "MissingData"


 #Generates a high resolution reference space to sample thickness
 c3d $WDIR/${what}_seg.nii.gz -trim 3x3x3vox -interpolation NearestNeighbor \
 -resample 300% -scale 0 -type uchar -o $WDIR/ref.mha
 
 #Fill the tetmesh. by default, background is zero. Do this when normalizing by reference
 $C3D_HOME/tetfill -c "VoronoiRadius" $SUBJ_TETRA $WDIR/ref.mha $WDIR/${id}_${what}_fillthk.mha

 ## Dilate the thickness image by a few pixels so mesh_image_sample is less sensitive to outlier vertices
 c3d $WDIR/${id}_${what}_fillthk.mha -thresh 0.001 inf 1 0 -o $TMPDIR/${id}_${what}_mask.mha

# To check image space compatability
# $C3D_HOME/c3d -verbose $SUBJSEG $WDIR/${what}_seg.nii.gz -times -trim 3x3x3vox -interpolation NearestNeighbor -resample 300% -info \
# $WDIR/${id}_${what}_fillthk.mha -info

 # Take the segmentation, retain only the labels pertaining to the list of interest
 # Use this to assign missing labels to the thickness image
 $C3D_HOME/c3d $SUBJSEG $WDIR/${what}_seg.nii.gz -times -trim 3x3x3vox -interpolation NearestNeighbor -resample 300% \
 -retain-labels 1 2 -replace 2 NaN $WDIR/${id}_${what}_fillthk.mha -times -o $THKIMG_MASKED 
#  \
#  -thresh 0.001 inf 1 0 -o $WDIR/${id}_${what}_thkmask_md.mha 

<<"SKIPOLD"
 # Sample the mesh by the segmentation
 $ROOT/pkgs/mesh_image_sample -B -b $RESLICED_SKEL $THKIMG_MASKED $THK_SKEL Thickness_nan

 $ROOT/pkgs/mesh_image_sample -B -b $THK_SKEL $WDIR/${id}_${what}_thkmask_md.mha $THK_SKEL Thickness_ref

 python NormalizeSkelThickness.py -i $THK_SKEL -o $THK_SKEL -d "Thickness_nan" -r "Thickness_ref"
SKIPOLD


  if [[ $what == "gm" ]]; then

     #Also add the histology label to the thickness mesh here!
    SUBJSEG=$WORK/template_tosubj/${id}/${id}_unwarp_histoseg_coarse.nii.gz

    c3d $SUBJSEG -replace 5 0 13 0 -o $TMPDIR/histo_seg_coarse_${what}.nii.gz

    # Take the segmentation, retain only the labels pertaining to the list of interest
    # expand these labels and use them to assign a label to each vertex on the thickness
    # mesh.
    c3d -verbose $TMPDIR/histo_seg_coarse_${what}.nii.gz -trim 5vox -replace 0 1000 \
      -split \
      -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
      -scale 0 -merge \
      -replace 1000 0 -o $TMPDIR/seg_expanded.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $THK_SKEL $TMPDIR/seg_expanded.nii.gz $WDIR/${id}_${what}_warped_thickness_skeleton_whistolabels.vtk histo_label_coarse

    #Fine histology annotation
    SUBJSEG=$WORK/template_tosubj/${id}/${id}_unwarp_histoseg_fine.nii.gz

    c3d $SUBJSEG -replace 5 0 -o $TMPDIR/histo_seg_fine_${what}.nii.gz

    # Take the segmentation, retain only the labels pertaining to the list of interest
    # expand these labels and use them to assign a label to each vertex on the thickness
    # mesh.
    c3d -verbose $TMPDIR/histo_seg_fine_${what}.nii.gz -trim 5vox -replace 0 1000 \
      -split \
      -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
      -scale 0 -merge \
      -replace 1000 0 -o $TMPDIR/seg_expanded.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $WDIR/${id}_${what}_warped_thickness_skeleton_whistolabels.vtk $TMPDIR/seg_expanded.nii.gz $WDIR/${id}_${what}_warped_thickness_skeleton_whistolabels.vtk histo_label_fine

  fi


<<"FIN"
  #For image space GLM, map each subjects' thickness image into template space
  #Also map the mask
  greedy -d 3 \
  -rf $SEG \
  -rm $THKIMG_MASKED $WDIR/${id}_${what}_thickness_md_reslicedtotemplate.nii.gz \
  -ri LABEL 0.24vox \
  -rm $WDIR/${id}_${what}_thkmask_md.mha $WDIR/${id}_${what}_thickness_md_binmask_reslicedtotemplate.nii.gz \
  -r $(cat $SDIR/chain_unwarped_to_template.txt)

  #Apply smoothing
  c3d $WDIR/${id}_${what}_thickness_md_reslicedtotemplate.nii.gz -smooth-fast 0.6mm -as T \
  $WDIR/${id}_${what}_thickness_md_binmask_reslicedtotemplate.nii.gz -as M -smooth-fast 0.6mm -shift 0.00001 -reciprocal \
  -push M -times -push T -times -o $WDIR/${id}_${what}_thickness_md_smoothed06_reslicedtotemplate.nii.gz

FIN

}


# Correlation with tau (excluding tdp subjects) but covarying ONLY for age
function thickness_glm_tau()
{

  # This function performs GLM on dark band thickness
  what=${1?}

  # The group membership file
  GRPFILE=$ROOT/scripts/path_data/groups_taucorr_contra_v3_2022_NFTquant.txt
  #GRPFILE=$ROOT/scripts/path_data/groups_taucorr_contra_v3_2022.txt

  # The template thickness mesh - SHOULD USE TSM
  REFMESH=$THKDIR/template_${what}_tsm_skel.vtk

  # Work directory for this analysis
  WDIR=$THKDIR/glm_tau/${what}_thickness_age_sex

  # The input mesh for GLM
  GLMINPUT_MD=$WDIR/glm_input_md_skel.vtk
  GLMRESULT=$WDIR/glm_result_skel.vtk
  GLMOUTPUT=$WDIR/glmoutput_skel.txt
  
  # Create the directory
  mkdir -p $WDIR

  # List of meshes to average
  local MERGE_LIST=""

  # Clear the design matrix
  DESIGN_MAT=$WDIR/design.txt
  rm -rf $DESIGN_MAT

  for subj in $(tail -n +2 $GRPFILE | grep -v 'NO' | awk '{print $1}'); do

    # Get the design matrix info
    local stage age dummy tdp ab asyn A B C sex neuronloss q_tangles q_threads NFT_90_35 NFT_MEAN_35 NFT_90_CA NFT_MEAN_CA NFT_90_EC NFT_MEAN_EC
    read dummy stage tdp asyn ab age sex A B C neuronloss q_tangles q_threads NFT_90_35 NFT_MEAN_35 NFT_90_CA NFT_MEAN_CA NFT_90_EC NFT_MEAN_EC <<< $(cat $GRPFILE | grep $subj)
    #if [[ $age -ge 50 ]] && [[ $(echo $NFT_90_35 | sed 's/\r//') != 'NA' ]]; then #&& [[ $subj != 'HNL29_18-L' ]]
    if [[ $age -ge 50 ]]; then #&& [[ $subj != 'HNL29_18-L' ]]
	 

      if [[ $sex == "M" ]]; then
	sex_var=0
      else
	sex_var=1
      fi

      echo $subj

      local MERGE_ENTRY=$(ls $THKDIR/subj_skel/${subj}/*${what}_warped_thickness_skeleton.vtk)

      # Generate the design matrix value
      echo $stage $age $sex_var 1 >> $DESIGN_MAT

      # Append to the merge list
      MERGE_LIST="$MERGE_LIST $MERGE_ENTRY"
    fi

  done

  # Combine meshes into a single VTK mesh for analysis
  $ROOT/pkgs/cmrep_vtk9/mesh_merge_arrays -B -r $REFMESH $GLMINPUT_MD Thickness_nan $MERGE_LIST

   # Generate the contrast and nuissance
  echo "-1 0 0 0" > $WDIR/contrast.txt
  echo "0 1 1 0" > $WDIR/nuissance.txt

 # Perform GLM with permutation - freedman-lane 
 # Added -T option to triangulate mesh  
  # $ROOT/pkgs/meshglm -m $GLMINPUT_MD $GLMRESULT -a Thickness_nan -d 4 \
  # -g $DESIGN_MAT $WDIR/contrast.txt -f -n $WDIR/nuissance.txt \
  # -T -M 5 -p 1000 -s P -t 0.01 -e > $GLMOUTPUT

  #Now specifying that a min of 5 observations at each vertex should not be NaN
  $ROOT/pkgs/meshglm -m $GLMINPUT_MD $GLMRESULT -a Thickness_nan -d 4 \
  -g $DESIGN_MAT $WDIR/contrast.txt \
  -T -M 5 -p 1000 -s P -t 0.01 -e > $GLMOUTPUT

 if [[ $what == "gm" ]]; then

  if [[ ! -f $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz ]]; then
    #Probability of each voxel's label is based on its distance to that region. -dup -times -sqrt computes the absolute value of the distance
    c3d $ADIR/histology_annotation_coarse/template_hf_label.nii.gz -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz
  fi

    $ROOT/pkgs/mesh_image_sample -B -i 0 $GLMRESULT $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz $WDIR/glm_result_skel_hf.vtk HF_label

    #Full histology annotation
    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -replace 5 0 13 0 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $WDIR/glm_result_skel_hf.vtk $TMPDIR/seg_expanded_coarse.nii.gz $WDIR/glm_result_skel_hf.vtk histo_label_coarse


  elif [[ $what == "srlm" ]]; then

    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -retain-labels 5 13 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $GLMRESULT $TMPDIR/seg_expanded_coarse.nii.gz $WDIR/glm_result_skel_hf.vtk histo_label_coarse


 fi


}

#Correlation with tau but covarying for TDP43
function thickness_glm_tau_tdp()
{
  # This function performs GLM on gm and dark band thickness
  what=${1?}

  # The group membership file
  GRPFILE=$ROOT/scripts/path_data/groups_taucorr_contra_v3_2022.txt

  # The template thickness mesh - SHOULD USE TSM
 REFMESH=$THKDIR/template_${what}_tsm_skel.vtk

  # Work directory for this analysis
  WDIR=$THKDIR/glm_asyn/${what}_thickness_tdp_tau_age_sex

  # The input mesh for GLM
  GLMINPUT_MD=$WDIR/glm_input_md_skel.vtk
  GLMRESULT=$WDIR/glm_resul_skel.vtk
  GLMOUTPUT=$WDIR/glmoutput_skel.txt
  
  # Create the directory
  mkdir -p $WDIR

  # List of meshes to average
  local MERGE_LIST=""

  # Clear the design matrix
  DESIGN_MAT=$WDIR/design.txt

  rm -rf $DESIGN_MAT

  for subj in $(tail -n +2 $GRPFILE | grep -v 'OTHER' | awk '{print $1}'); do

    # Get the design matrix info
    local stage age dummy tdp ab asyn A B C sex neuronloss sustage braak q_tangle q_threads
    read dummy stage tdp asyn ab age sex A B C neuronloss q_tangle q_threads <<< $(cat $GRPFILE | grep $subj)

    if [[ $age -ge 50 ]]; then
   
      local MERGE_ENTRY=$(ls $THKDIR/subj_skel/${subj}/*${what}_warped_thickness_skeleton.vtk)

	if [[ $sex == "M" ]]; then
		sex_var=0
	else
		sex_var=1
	fi

    	# Generate the design matrix containing tdp value
    	echo $stage $tdp $asyn $age $sex_var 1 >> $DESIGN_MAT

    	# Append to the merge list
    	MERGE_LIST="$MERGE_LIST $MERGE_ENTRY"
    fi

  done

  # Combine meshes into a single VTK mesh for analysis - was "Thickness_nan_norm" before
  $ROOT/pkgs/cmrep_vtk9/mesh_merge_arrays -B -r $REFMESH $GLMINPUT_MD Thickness_nan $MERGE_LIST

  # #Generate the contrast and nuissance containing tdp
  echo "0 0 -1 0 0 0" > $WDIR/contrast.txt
  echo "1 1 0 1 1 0" > $WDIR/nuissance.txt

  # Perform GLM with permutation - freedman-lane and co-vary for tdp
  # Original diffusion was 4, but reduced to 1 for the newer results for neater cluster.
  $ROOT/pkgs/meshglm -m $GLMINPUT_MD $GLMRESULT -a Thickness_nan -d 4 \
  -g $DESIGN_MAT $WDIR/contrast.txt -f -n $WDIR/nuissance.txt \
  -p 1000 -M 0.5 -s P -t 0.01 -T -e --threads 4 > $GLMOUTPUT

 if [[ $what == "gm" ]]; then

    $ROOT/pkgs/mesh_image_sample -B -i 0 $GLMRESULT $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz $WDIR/glm_result_skel_hf.vtk HF_label

    #Full histology annotation
    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -replace 5 0 13 0 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $WDIR/glm_result_skel_hf.vtk $TMPDIR/seg_expanded_coarse.nii.gz $WDIR/glm_result_skel_hf.vtk histo_label_coarse


  elif [[ $what == "srlm" ]]; then

    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -retain-labels 5 13 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $GLMRESULT $TMPDIR/seg_expanded_coarse.nii.gz $WDIR/glm_result_skel_hf.vtk histo_label_coarse


 fi

}

function thickness_glm_age()
{
  # This function performs GLM on dark band thickness
  what=${1?}

  # The group membership file
  GRPFILE=$ROOT/scripts/path_data/groups_taucorr_contra_v3_2022.txt

  # The template thickness mesh - SHOULD USE TSM
  REFMESH=$THKDIR/template_${what}_tsm_skel.vtk

  # Work directory for this analysis
  WDIR=$THKDIR/glm_age/${what}_thickness

  # The input mesh for GLM
  GLMINPUT_MD=$WDIR/glm_input_md_skel.vtk
  GLMRESULT=$WDIR/glm_result_skel.vtk
  GLMOUTPUT=$WDIR/glmoutput_skel.txt

  # Create the directory
  mkdir -p $WDIR

  # List of meshes to average
  local MERGE_LIST=""

  # Clear the design matrix
  rm -rf $WDIR/design.txt

  # Iterate over all subjects in the group file
  for subj in $(tail -n +2 $GRPFILE | grep -v 'OTHER' | awk '{print $1}'); do

    # Get the design matrix info
    local stage age dummy tdp ab asyn A B C sex neuronloss sustage braak q_tangle q_threads
    read dummy stage tdp asyn ab age sex A B C neuronloss q_tangle q_threads <<< $(cat $GRPFILE | grep $subj)

    if [[ $age -ge 50 ]]; then

      local MERGE_ENTRY=$(ls $THKDIR/subj_skel/${subj}/*${what}_warped_thickness_skeleton.vtk)

    	# Append to the merge list
    	MERGE_LIST="$MERGE_LIST $MERGE_ENTRY"

    	# Generate the design matrix value
    	echo 1 $age >> $WDIR/design.txt
    fi
    
  done

  # Combine meshes into a single VTK mesh for analysis
  $ROOT/pkgs/cmrep_vtk9/mesh_merge_arrays -B -r $REFMESH $GLMINPUT_MD Thickness_nan $MERGE_LIST

  # Generate the contrast 
  echo "0 -1" > $WDIR/contrast.txt

  # Perform GLM with permutation - Need the -T to perform triangulation. Clears the mesh of any vtkLines
  $ROOT/pkgs/meshglm -m $GLMINPUT_MD $GLMRESULT -a Thickness_nan -d 4 -T \
  -g $WDIR/design.txt $WDIR/contrast.txt -M 5 -p 1000 -s P -t 0.01 --threads 4 -e > $GLMOUTPUT

  
 if [[ $what == "gm" ]]; then

    if [[ ! -f $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz ]]; then
      #Probability of each voxel's label is based on its distance to that region. -dup -times -sqrt computes the absolute value of the distance
      c3d $ADIR/histology_annotation_coarse/template_hf_label.nii.gz -trim 5vox -replace 0 1000 \
      -split \
      -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
      -scale 0 -merge \
      -replace 1000 0 -o $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz
    fi

    $ROOT/pkgs/mesh_image_sample -B -i 0 $GLMRESULT $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz $WDIR/glm_result_skel_hf.vtk HF_label

    #Full histology annotation
    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -replace 5 0 13 0 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $WDIR/glm_result_skel_hf.vtk $TMPDIR/seg_expanded_coarse.nii.gz $WDIR/glm_result_skel_hf.vtk histo_label_coarse


  elif [[ $what == "srlm" ]]; then

    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -retain-labels 5 13 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $GLMRESULT $TMPDIR/seg_expanded_coarse.nii.gz $WDIR/glm_result_skel_hf.vtk histo_label_coarse


 fi

}

function thickness_glm_misc()
{
  # This function performs GLM on gm andn dark band thickness
  what=${1?}

  independent=${2?}

  # The group membership file
  #GRPFILE=$ROOT/scripts/path_data/groups_taucorr_contra_v3_2022.txt
  GRPFILE=$ROOT/scripts/path_data/groups_taucorr_contra_v3_2022_NFTquant.txt

  # The template thickness mesh - SHOULD USE TSM
  REFMESH=$THKDIR/template_${what}_tsm_skel.vtk

  # Work directory for this analysis
  WDIR=$THKDIR/glm_QuantTau/${what}_thickness_${independent}_age_sex

  # The input mesh for GLM
  GLMRESULT=$WDIR/glm_result_skel.vtk
  GLMINPUT_MD=$WDIR/glm_input_md_skel.vtk
  GLMOUTPUT=$WDIR/glmoutput_skel.txt

  # Create the directory
  mkdir -p $WDIR

  # List of meshes to average
  local MERGE_LIST=""
  
  # Clear the design matrix
  rm -rf $WDIR/design.txt

  # Iterate over all subjects in the group file
  for subj in $(tail -n +2 $GRPFILE | grep -v 'OTHER' | awk '{print $1}'); do

    # Get the design matrix info
    local stage age dummy tdp ab asyn A B C sex neuronloss q_tangles q_threads NFT_MEAN_35 NFT_MEAN_CA NFT_MEAN_EC NFT_MEAN_PrPaS NFT_MEAN_TE
    read dummy stage tdp asyn ab age sex A B C neuronloss q_tangles q_threads NFT_MEAN_35 NFT_MEAN_CA NFT_MEAN_EC NFT_MEAN_PrPaS NFT_MEAN_TE <<< $(cat $GRPFILE | grep $subj)

    if [[ $independent == 'neuronloss' ]]; then
      var_interest=$neuronloss
    elif [[ $independent == 'sustage' ]]; then
    	var_interest=$sustage
    elif [[ $independent == 'braak' ]]; then
    	var_interest=$braak
    elif [[ $independent == 'q_thread' ]]; then
      var_interest=$q_threads
    elif [[ $independent == 'q_tangle' ]]; then
      var_interest=$q_tangles
    elif [[ $independent == 'NFT_MEAN_EC' ]]; then
      var_interest=$NFT_MEAN_EC
    elif [[ $independent == 'NFT_MEAN_BA35' ]]; then
      var_interest=$NFT_MEAN_35
    elif [[ $independent == 'NFT_MEAN_CA1' ]]; then
      var_interest=$NFT_MEAN_CA
    elif [[ $independent == 'NFT_MEAN_PrPaS' ]]; then
      var_interest=$NFT_MEAN_PrPaS
    elif [[ $independent == 'NFT_MEAN_TE' ]]; then
      var_interest=$NFT_MEAN_TE
    fi

    echo $var_interest
    if [[ $(echo $var_interest | sed 's/\r//') != 'NA' ]]; then

      if [[ $sex == "M" ]]; then
	sex_var=0
      else
	sex_var=1
      fi

      local MERGE_ENTRY=$(ls $THKDIR/subj_skel/${subj}/*${what}_warped_thickness_skeleton.vtk)
    
      # Append to the merge list
      MERGE_LIST="$MERGE_LIST $MERGE_ENTRY"

      # Generate the design matrix value
      echo 1 $var_interest $age $sex_var >> $WDIR/design.txt
    fi
  done

 # Combine meshes into a single VTK mesh for analysis - was "Thickness_nan_norm" before
  $ROOT/pkgs/cmrep_vtk9/mesh_merge_arrays -B -r $REFMESH $GLMINPUT_MD Thickness_nan $MERGE_LIST

  # #Generate the contrast and nuissance containing tdp
  echo "0 -1 0 0" > $WDIR/contrast.txt
  echo "0 0 1 1" > $WDIR/nuissance.txt

  #Remove -M 5 flag for high/low Braak analysis - was getting segmentation fault 
  # Perform GLM with permutation - freedman-lane and co-vary for tdp
  # Original diffusion was 4, but reduced to 1 for the newer results for neater cluster.
  $ROOT/pkgs/meshglm -m $GLMINPUT_MD $GLMRESULT -a Thickness_nan -d 4 -T \
  -g $WDIR/design.txt $WDIR/contrast.txt -f -n $WDIR/nuissance.txt -M 0.5 -p 1000 -s P -t 0.01 --threads 4 -e > $GLMOUTPUT

 if [[ $what == "gm" ]]; then

    $ROOT/pkgs/mesh_image_sample -B -i 0 $GLMRESULT $ADIR/histology_annotation_coarse/seg_expanded_hflabel.nii.gz $WDIR/glm_result_skel_hf.vtk HF_label

    #Full histology annotation
    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -replace 5 0 13 0 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $WDIR/glm_result_skel_hf.vtk $TMPDIR/seg_expanded_coarse.nii.gz $WDIR/glm_result_skel_hf.vtk histo_label_coarse


  elif [[ $what == "srlm" ]]; then

    c3d -verbose $ADIR/histology_annotation_coarse/template_histo_seg_final.nii.gz -retain-labels 5 13 -trim 5vox -replace 0 1000 \
    -split \
    -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    -scale 0 -merge \
    -replace 1000 0 -o $TMPDIR/seg_expanded_coarse.nii.gz

    # Sample the mesh by the segmentation
    $ROOT/pkgs/mesh_image_sample -B -i 0 $GLMRESULT $TMPDIR/seg_expanded_coarse.nii.gz $WDIR/glm_result_skel_hf.vtk histo_label_coarse


 fi

}

function thickness_hotspot_summary()
{

for id in $(cat $SUBJ_LIST); do

    local WDIR=$WORK/thickness/subj/${id}

    cat $WDIR/${id}_gm_thickness_by_hotspot.txt \
     $WDIR/${id}_srlm_thickness_by_hotspot.txt \
      > $WDIR/${id}_all_thickness_by_hotspot.txt
    sed -i '1d;4d' $WDIR/${id}_all_thickness_by_hotspot.txt
    sed -i 's/\s\+/,/g' $WDIR/${id}_all_thickness_by_hotspot.txt
 done

 TAU_90_CSV=$ADIR/tau_atlas/template_density_Tau_tangles_hotspots_90_tiles.csv
 for id in $(cat $ADIR/scripts/subj_tau_atlas.txt | awk '{print $1}'); do
      local WDIR=$WORK/thickness/subj/${id}
      cat $TAU_90_CSV | grep $id > $WDIR/${id}_all_tau_by_hotspot.txt
      paste $WDIR/${id}_all_thickness_by_hotspot.txt $WDIR/${id}_all_tau_by_hotspot.txt
   done > $WORK/thickness/thickness_tau_by_hotspot.txt

   for id in $(cat $SUBJ_LIST | awk '{print $1}'); do
      local WDIR=$WORK/thickness/subj/${id}
      cat $WDIR/${id}_all_thickness_by_hotspot.txt
   done > $WORK/thickness/thickness_by_hotspot_allsubj.txt

}

function thickness_ROI_summary()
{

  WDIR="$THKDIR/thickness_ROI"

  for id in $(cat $SUBJ_LIST); do

    cat $WDIR/${id}_srlm_lstat_thickness_coarse.txt \
        | awk -v id=$id 'NR > 2 {print id,$1,$2}' 
  done > $WDIR/all_thickness_by_srlm_histolabel.txt

  for id in $(cat $SUBJ_LIST); do

    WDIR=$THKDIR/thickness_ROI

    cat $WDIR/${id}_gm_lstat_thickness_coarse.txt \
        | awk -v id=$id 'NR > 2 {print id,$1,$2}' 
  done > $WDIR/all_thickness_by_gm_histolabel.txt

 cat \
      $WDIR/all_thickness_by_srlm_histolabel.txt \
      $WDIR/all_thickness_by_gm_histolabel.txt \
      > $WDIR/all_thickness_by_histolabel.txt

for id in $(cat $SUBJ_LIST); do

    cat $WDIR/${id}_gm_lstat_thickness_ECsubfields.txt \
        | awk -v id=$id 'NR > 2 {print id,$1,$2}' 
  done > $WDIR/all_thickness_by_ECsubfields_histolabel.txt

<<"SKIP"
    cat \
      $WDIR/${id}_gm_thickness_by_label.txt \
      $WDIR/${id}_srlm_thickness_by_label.txt \
      > $WDIR/${id}_all_thickness_by_label.txt

    # Note the multiplication by two here
    echo ${id} $(for label in 1 2; \
      do cat $WDIR/${id}_all_thickness_by_label.txt \
        | grep "^${id} ${label} " | awk '{print $3 * 2.0}'; done) \
        > $WDIR/${id}_all_thickness_by_label_wide.txt

    cat \
      $WDIR/${id}_gm_thickness_by_histolabel.txt \
      $WDIR/${id}_srlm_thickness_by_histolabel.txt \
      > $WDIR/${id}_all_thickness_by_histolabel.txt 
 
 done

   for id in $(cat $SUBJ_LIST); do

      if [[ $mode == "template" ]]; then
        WDIR=$THKDIR/subj/${id}
       elif [[ $mode == "native" ]]; then
        WDIR=$THKDIR/subj_thickness_unwarp/${id}
       fi

      cat $WDIR/${id}_all_thickness_by_histolabel.txt 
   done > $THKDIR/thickness_by_histolabel_${mode}.txt

SKIP

}

function map_seg_to_subject()
{

  id=${1?}

  mode=${2?}

  # Map the atlas segmentation back into subject space
  WDIR=$WORK/template_tosubj/${id}
  SDIR=$ADIR/mst_multilabel/final/${id}

  # Unwarped image
  IMGUW=$SDIR/${id}_unwarp_img.nii.gz 

  IMG_AXISALIGN=/project/hippogang_3/sravikumar/atlasPHG2019/pmatlas_raw/${id}/hires_n4_axisalign/${id}_n4clip_axisalign.nii.gz
  
 if [[ ! -f $IMG_AXISALIGN ]]; then
  	IMG_AXISALIGN=/project/hippogang_3/sravikumar/atlasPHG2019/preproc/${id}/${id}_axisalign_img.nii.gz
 fi

  HIRES_MRI_MANUAL_PHGSEG_AFFINE=/project/hippogang_3/sravikumar/atlasPHG2019/preproc/${id}/${id}_transform_to_axisalign.mat
  UNWARP_MAT=/project/hippogang_3/sravikumar/atlasPHG2019/preproc/${id}/unwarp/${id}_unwarp.mat
   
  # The segmentation we are going to use
  SEG=$ADIR/histology_annotation_${mode}/template_histo_seg_final.nii.gz
  
  mkdir -p $WDIR

  # if [[ ! -f $(ls $WDIR/${id}_unwarp_img*.nii.gz) ]]; then
  #     if [[ $id == "HNL44_19-L" ]]; then
  #       #Pad the reference unwarped image since it cuts off the MTL
  #       c3d $IMGUW -pad 30x30x30vox 30x30x30vox -o $WDIR/${id}_unwarp_img_padded.nii.gz
  #     else
  #       ln -sf $IMGUW $WDIR/
  #     fi
  # fi

  # IMG_REF=$(ls $WDIR/${id}_unwarp_img.nii.gz)

<<"DONE" 
IMG_REF=$THKDIR/subj/$id/${id}_unwarp_img_padded.nii.gz

$GREEDY_HOME/greedy -d 3 \
  -ri LABEL 0.24vox \
  -rf $IMG_REF  \
  -rm $SEG $WDIR/${id}_unwarp_histoseg_${mode}.nii.gz \
  -r $(cat $SDIR/chain_template_to_unwarped.txt)
DONE

$GREEDY_HOME/greedy -d 3 \
  -ri LABEL 0.24vox \
  -rf $IMG_AXISALIGN  \
  -rm $SEG $WDIR/${id}_axisalign_histoseg_${mode}.nii.gz \
  -r $HIRES_MRI_MANUAL_PHGSEG_AFFINE $UNWARP_MAT,-1 $(cat $SDIR/chain_template_to_unwarped.txt)
 
# c3d $IMG_REF $WDIR/${id}_unwarp_histoseg_${mode}.nii.gz \
#    -lstat > $WDIR/${id}_lstat_${mode}.txt

}

function mask_badreg_thickness()
{

	#VTK library needed for mesh_fill_missing_data
	export LD_LIBRARY_PATH=/project/hippogang_3/sravikumar/packages/vtk-release-6.3/bin/lib

	id=${1?}

	WDIR=$THKDIR/subj/$id

	#Registration mask is based on final geodesic shooting warp
	GSHOOT_WARP=$ADIR/mst_multilabel/gshoot/$id/iter_4/shooting_warp.nii.gz
	
	#Reference unwarp space for subject
	SDIR=$ADIR/mst_multilabel/final/${id}
  IMGUW=$SDIR/${id}_unwarp_img.nii.gz

	#Subject segmentation with missing label in subject space
  #April 2023 - the unwarped image cuts of part of the segmentation. Instead of changing the original 
  # unwarped image (might have transformations?), just updated the segmentation here that is used for thickness 
  # measurements
  #April 2023 - do this for all cases not just 120126 (120126 was only padded with 30x0x0vox 0x0x0vox)
  if [[ ! -f $WDIR/${id}_unwarp_img_padded.nii.gz ]]; then
    c3d $IMGUW -pad 30x30x30vox 30x30x30vox -o $WDIR/${id}_unwarp_img_padded.nii.gz
  fi

  IMGUW=$WDIR/${id}_unwarp_img_padded.nii.gz

  greedy -d 3 \
  -rf $IMGUW \
  -ri LABEL 0.24vox \
  -rm $ROOT/preproc/$id/${id}_axisalign_phgseg_singlelabel.nii.gz $WDIR/${id}_unwarp_phgseg_singlelabel_padded.nii.gz \
  -r $ROOT/preproc/$id/unwarp/${id}_unwarp.mat $ROOT/preproc/$id/${id}_transform_to_axisalign.mat,-1

  SUBJSEG=$WDIR/${id}_unwarp_phgseg_singlelabel_padded.nii.gz

	#Output bad-registration masks
	GSHOOT_WARP_MASK=$WDIR/${id}_shooting_999_warpmask.nii.gz
	GSHOOT_WARP_MASK_UNWARP=$WDIR/${id}_unwarp_shooting_999_warpmask.nii.gz
	MISSING_LABEL_MASK_UNWARP=$WDIR/${id}_unwarp_missinglabel_withbadreg.nii.gz	

	#Thickness meshes
	local THKMESH_MD=$WDIR/${id}_gm_thickness_md_badreg.vtk
  local THKMESH=$WDIR/${id}_gm_thickness.vtk

	#if [[ ! -f $GSHOOT_WARP_MASK_UNWARP ]]; then

		#Create a mask of regions with deformation magnitude in the top 0.1%
		c3d -pim fq -mcs $GSHOOT_WARP \
		-foreach -info -dup -times -endfor \
		-accum -add -endaccum -sqrt \
    -thresh 99.9% inf 1 0 -o $GSHOOT_WARP_MASK

		#Warp to subject space, and combine with artifact label. To specify NaNs
		$GREEDY_HOME/greedy -d 3 \
    -ri LABEL 0.24vox \
    -rf $IMGUW \
    -rm $GSHOOT_WARP_MASK $GSHOOT_WARP_MASK_UNWARP \
    -r $(cat $SDIR/chain_template_to_unwarped.txt)
	#fi

	#Combine registration mask with artifact label
	c3d $SUBJSEG -as M -thresh 1 inf 1 0 \
	$GSHOOT_WARP_MASK_UNWARP -times \
	-push M -add -replace 3 2 \
	-o $MISSING_LABEL_MASK_UNWARP

	# Expand the segmentation
  c3d $MISSING_LABEL_MASK_UNWARP -retain-labels 1 2 -trim 5vox -replace 0 1000 \
  -split \
  -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
 	-scale 0 -merge \
 	-replace 1000 0 -o $TMPDIR/seg_expanded.nii.gz

  # Sample the mesh by the segmentation
  $CMREP_HOME/mesh_image_sample -i 0 $THKMESH $TMPDIR/seg_expanded.nii.gz $THKMESH LABEL_BADREG

  # Use SR script to replace missing data points with NaN. Outputs vtk mesh in binary format
  # Do this after I sample the histo ROI instead -since mesh image sample doesn't read in binary vtk file
 	$MESH_SRCODE/mesh_fill_missing_data -b $THKMESH Thickness LABEL_BADREG 2 $THKMESH_MD

}

function thickness_ROI_stats()
{

<<"FIN"
N=$(cat $SUBJ_LIST | wc -l)

# Generate csv/text files with thickness measures per ROI
for what in gm srlm; do 
 	for ((i=1;i<=$N;i++)); do

 	  subj=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
 	       
    #Use histology labels to generate ROI mean/median thickness
    $ROOT/scripts/job_launcher.sh -m 25G -o $WORK/dump -N "thk_label_${what}_${subj}" \
    $0 thickness_perlabel $subj $what 

 	done
done

# Generate csv/text files with thickness measures per ROI
for ((i=1;i<=$N;i++)); do

 	subj=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
 	       
  #Use histology labels to generate ROI mean/median thickness
  $ROOT/scripts/job_launcher.sh -m 25G -o $WORK/dump -N "thk_EClabel_${subj}" \
  $0 thickness_perlabel_EC $subj 

done

$ROOT/scripts/job_launcher.sh -w "thk_EClabel_*" /bin/sleep 0

FIN

thickness_ROI_summary 

}


function meshglm_analyses(){

#Perform thickness GLM - Modified Nov 2022 to perform glm on skeleton
for what in gm srlm; do #gm
	
  #Thickness vs age
  #$ROOT/scripts/job_launcher.sh -m 30G -o $WORK/dump -N "glm_age_${what}" $ROOT/scripts/template_skel_analysis.sh thickness_glm_age ${what}

  # Thickness vs Tau - Age is the only covariate. Run this using the NFT only subset as well
  #$ROOT/scripts/job_launcher.sh -m 60G -o $WORK/dump -N "glm_corr_tau_${what}" $ROOT/scripts/template_skel_analysis.sh thickness_glm_tau ${what}

  #Correct for TDP-43 and aSyn and age. Also run with TDP-43 and asyn as the variable of interest
  $ROOT/scripts/job_launcher.sh -m 60G -o $WORK/dump -N "glm_corr_tau_copath_${what}" $ROOT/scripts/template_skel_analysis.sh thickness_glm_tau_tdp ${what} 

  for mode in NFT_MEAN_BA35 NFT_MEAN_EC NFT_MEAN_CA1 NFT_MEAN_PrPaS NFT_MEAN_TE; do
		  
      $ROOT/scripts/job_launcher.sh -m 30G -o $WORK/dump -N \
	"glm_corr_${what}_${mode}" $ROOT/scripts/template_skel_analysis.sh \
	thickness_glm_misc ${what} $mode
   done 

done

}

function main()
{
  
  mkdir -p $WORK/dump

  N=$(cat $SUBJ_LIST | wc -l)

  # THICKNESS ANALYSIS
  # extract meshes from template (GM and SRLM)
  #$ROOT/scripts/job_launcher.sh -m 15G -o $WORK/dump -N "thk_templates" "$0" thickness_templates

  # $ROOT/scripts/job_launcher.sh -w "thk_templates" /bin/sleep 0
  
<<"DONE"
  #Create a segmentation mask in each subject's space masking out regions of poor registration quality 
  #- already applied to skel in previous step. Created in old script
  #Rerun this for 120126 (N - 32)- pad unwarped image since part of ERC is being cut off  (April 2023)
  # Parts of HNNL32, 35 and 44 are all getting cut off so just pad all cases to be safe
  for ((i=1;i<=$N;i++)); do
    subj=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
  	$ROOT/scripts/job_launcher.sh -m 15G -o $WORK/dump -N "mask_badreg_${subj}" $ROOT/scripts/template_skel_analysis.sh mask_badreg_thickness $subj
  done

  # Map the atlas histology segmentation into the individual subject spaces
  for subj in $(cat $SUBJ_LIST); do

  $ROOT/scripts/job_launcher.sh -m 20G -o $WORK/dump -N "map_coarse_histoseg_to_${subj}" \
  "$0" map_seg_to_subject $subj "coarse"

#  $ROOT/scripts/job_launcher.sh -m 30G -o $WORK/dump -N "map_fine_histoseg_to_${subj}" \
#  "$0" map_seg_to_subject $subj "fine"

  done
DONE

 # $ROOT/scripts/job_launcher.sh -w "map_*_histoseg_to_*" /bin/sleep 0

<<"DONE"
  # compute thickness in each subject's space
  for what in gm srlm; do  #srlm
    for ((i=1;i<=$N;i++)); do
      
    subj=$(cat $SUBJ_LIST | head -n $i | tail -n 1)

    $ROOT/scripts/job_launcher.sh -m 50G -o $WORK/dump -N "thk_${what}_${subj}" \
    "$0" thickness_persubj_skel $subj $what
        
    done
  done

$ROOT/scripts/job_launcher.sh -w "thk_*" /bin/sleep 0
DONE

# Thickness statistical analysis - Need to do
#thickness_ROI_stats

#Need to run one more analysis just covarying for tdp/asyn and no age
meshglm_analyses


}


if [[ $1 ]]; then

  command=$1
  shift
	$command $@

else

  main

fi
