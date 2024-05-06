#!/bin/bash
#$ -S /bin/bash
set -x -e

ROOT=/project/hippogang_3/sravikumar/atlasPHG2019
source $ROOT/scripts/common.sh
#RAWDIR_OLD=/data/tesla-data/pauly/studies/pmatlas/ndri/import/nii
RAWDIR=/project/hippogang_3/sravikumar/atlasPHG2019/pmatlas_raw
ATLAS_2016=/project/hippogang_2/pauly/dan2015/atlas2016
CMREP_DIR=$ROOT/pkgs/cmrep-1.0.0-Linux-gcc64/bin/
#ANTS_PATH=/share/apps/ANTs/2014-06-23/build/bin/

# The directory where the preprocessing is done
WDIR=$ROOT/preproc

SUBJ_LIST=$ROOT/scripts/subj_atlas_2022_v2.txt

function organize_subj_2022()
{
  # This is the subject ID
  id=${1?}

  # Output for this work
  WORK=$WDIR/${id}
  mkdir -p $WORK/misc

  # The input image and manual segmentations - MTL and SRLM (as of 01/21/2021)
  INPUT_IMG=$ROOT/inputs/$id/${id}_axisalign_img.nii.gz

  #Naming is different for new cases
  if [[ ! -f $INPUT_IMG ]]; then
	INPUT_IMG=$ROOT/inputs/$id/${id}_n4clip_axisalign_img.nii.gz
  fi

  #Note: For new cases (incl. May 2022), the srlm segmentation is also in axisalign_phg (not a separate file) and
  # it includes the full extent of the SRLM - make sure all cleaned up segmentations have this filename
  INPUT_SEG=$ROOT/inputs/$id/${id}_axisalign_phg_fullsrlm.nii.gz

  #Corresponding ID according to PY naming scheme
  ID_PY=$(cat $ROOT/scripts/manual/hiresmri_manifest.txt \
    | awk -v id="$id" '$1 == id {print $2}')

  # Raw MRI scan - already n4 corrected and clipped
  #The n4 correction and clipping is now done when I copy files from flywheel (copy_inputs.sh).
  #So I don't need to do this step here anymore (March 2022)
  RAW_N4=$ROOT/pmatlas_raw/${ID_PY}/hires_n4/${id}_mri_raw_n4_clip.nii.gz
  N4_CLIP=$WORK/${id}_raw_n4clip.nii.gz
  
  ln -sf $RAW_N4 $N4_CLIP

 #Produce a link to the input data - i.e. the data segmented by Sadhana in the axisalign space
 # The segmentations contain multiple labels.
 #The green and yellow labels (2,4, 8) are used to mark regions which shouldn't be included in the intensity based registration.
 #Label 8 added here for SRLM artifacts. 6 is cysts in the SRLM. Don't include in intensity registration either
 #Can threshold to obtain single label segmentation for the shape matching step.
 #Specify filenames
  AXISALIGN_IMG=$WORK/${id}_axisalign_img.nii.gz
  AXISALIGN_SINGLE_WITHARTIFACTLBL=$WORK/${id}_axisalign_phgseg_singlelabel.nii.gz
  AXISALIGN_SINGLE=$WORK/${id}_axisalign_phgsegshape_singlelabel.nii.gz
  AXISALIGN_MULTI_WITHARTIFACTLBL=$WORK/${id}_axisalign_phgseg_multilabel.nii.gz
  AXISALIGN_MULTI=$WORK/${id}_axisalign_phgsegshape_multilabel.nii.gz
  AXISALIGN_SRLM=$WORK/${id}_axisalign_srlmseg_sr.nii.gz #Cysts and missing tissue are included as label 2
  AXISALIGN_SRLM_SHAPE=$WORK/${id}_axisalign_srlmsegshape_sr.nii.gz

  #Link to the MRI scan
  ln -sf $INPUT_IMG $AXISALIGN_IMG

  #Create version of the segmentation with and without the artifact labels.
  #Replace SRLM with GM for cases with combined input segmentaions
  $C3D_HOME/c3d $INPUT_SEG -replace 3 1 4 2 5 1 6 2 8 2 -o $AXISALIGN_SINGLE_WITHARTIFACTLBL
  $C3D_HOME/c3d $INPUT_SEG -thresh 1 8 1 0 -o $AXISALIGN_SINGLE
  $C3D_HOME/c3d $INPUT_SEG -replace 2 1 4 3 5 1 6 1 8 1 -o $AXISALIGN_MULTI

  #New cases don't have a separate srlm segmentation file - create a separate SRLM file
  #Inluded a check to make sure that the SRLM is a single connected component
  $C3D_HOME/c3d $INPUT_SEG -replace 5 1 6 2 8 2 -o $AXISALIGN_MULTI_WITHARTIFACTLBL
  
  $C3D_HOME/c3d $INPUT_SEG -retain-labels 5 6 8 -replace 5 1 6 2 8 2 -as SEG -comp -thresh 1 1 1 0 -o $AXISALIGN_SRLM_SHAPE -push SEG -times -o $AXISALIGN_SRLM

<<"SKIP"
 ## Create separate files for medial and lateral PHG seg
 $C3D_HOME/c3d $AXISALIGN_MULTI_WITHARTIFACTLBL -retain-labels 3 4 -replace 3 1 4 2 -o $WORK/${id}_axisalign_cslatseg.nii.gz
 $C3D_HOME/c3d $AXISALIGN_MULTI_WITHARTIFACTLBL -retain-labels 1 2 -o $WORK/${id}_axisalign_csmedseg.nii.gz
 $C3D_HOME/c3d $AXISALIGN_MULTI -thresh 3 3 1 0 -o $WORK/${id}_axisalign_cslatsegshape.nii.gz
 $C3D_HOME/c3d $AXISALIGN_MULTI -retain-labels 1 -o $WORK/${id}_axisalign_csmedsegshape.nii.gz
SKIP

UNWARP=$ROOT/phantom_unwarp_new/phantom_unwarp_affine_newparam.mat

# Create the unwarp dir
mkdir -p $WORK/unwarp
mkdir -p $ROOT/inputs/$id/segspace_transforms #This was needed for some of the old datasets, since I didn't manually generate the folder
if [[ ! -f $WORK/unwarp/${id}_unwarp.mat ]]; then
  cp -av $UNWARP $WORK/unwarp/${id}_unwarp.mat
fi

# Link to axisalign transformation
TRANSFORM_TO_AXISALIGN=$WORK/${id}_transform_to_axisalign.mat

# Manually generated transformation between raw and axisalign space. Segmentations are performed
# in axisalign space
MANUAL_RAW_TO_AXISALIGN=$ROOT/inputs/$id/segspace_transforms/${id}_raw_to_axisalign.mat
if [[ ! -f $MANUAL_RAW_TO_AXISALIGN ]]; then
	MANUAL_RAW_TO_AXISALIGN=$ROOT/scripts/manual/raw_to_axisalign_transforms/${id}_raw_to_axisalign.mat
fi

#For new subjects added in this dataset, segmentations were done in template space. Raw images were manually aligned with the template using ITK-SNAP registration tools:
# This matrix does not work for samples from the old dataset where I used the axisaligned images
if [[ "${id}" == "INDD119359-R" ]];then

        # Was oiginally using the wrong raw scan.Paul generated trasnformations
        # between the old and new raw scan
        ln -sf $MANUAL_RAW_TO_AXISALIGN  $TRANSFORM_TO_AXISALIGN

        REMAP_FOLDER=$ROOT/inputs/$id/INDD119359_remap
        CHAIN_SEG_TO_RAW="$REMAP_FOLDER/warp.nii.gz $REMAP_FOLDER/affine.mat $TRANSFORM_TO_AXISALIGN,-1"

elif [[ -f $MANUAL_RAW_TO_AXISALIGN ]]; then

   # Create a symbloic link to the raw to axis align transformation matrix
	 ln -sf $MANUAL_RAW_TO_AXISALIGN  $TRANSFORM_TO_AXISALIGN

 	 CHAIN_SEG_TO_RAW="$TRANSFORM_TO_AXISALIGN,-1"
 	 echo $id 2 >> $ROOT/scripts/seg_to_raw_atlas2022.txt

else

# If none of the above matrices exit,then this is a subject from the old atlas dataset.
# The transformation from unwrped raw_n4clip.nii.gz to axisalign space is stored in ../atlas2016/histolabels/axisalign/$id}/rigid/procrustes.mat

	# Produce a link to the tranformation from raw_n4clip to axisalign space
	cp --remove-destination $ATLAS_2016/histolabels/axisalign/$id/rigid/procrustes.mat \
 $WORK/${id}_transform_to_axisalign.mat

  CHAIN_SEG_TO_RAW="$WORK/${id}_warp.mat $TRANSFORM_TO_AXISALIGN,-1"
	echo $id 1 >> $ROOT/scripts/seg_to_raw_atlas2022.txt
fi

 #For HNL44, the raw image is cut off. Need to add before unwarping so that the missing tissue is included
 if [[ "${id}" == "HNL44_19-L" ]];then
	c3d $N4_CLIP -pad 0x0x0vox 0x0x20vox 0 -o /tmp/ref_space.nii.gz
	REF_SPACE=/tmp/ref_space.nii.gz
 else
	REF_SPACE=$N4_CLIP
 fi

 # Tranform the image and segmentation to the unwarped image space.
 # Use the clipped imageso that intensities are normalized across subjects
 # Intensity image
 $GREEDY_HOME/greedy -d 3 \
 -rf $REF_SPACE \
 -rm $N4_CLIP $WORK/unwarp/${id}_unwarp_img.nii.gz \
 -r $WORK/unwarp/${id}_unwarp.mat

 # Segmentation image - single label and multilabel - with and without artifact label
 $GREEDY_HOME/greedy -d 3 \
 -rf $WORK/unwarp/${id}_unwarp_img.nii.gz \
 -ri LABEL 0.24vox \
 -rm $AXISALIGN_SINGLE_WITHARTIFACTLBL $WORK/unwarp/${id}_unwarp_phgseg_singlelabel.nii.gz \
 -rm $AXISALIGN_MULTI_WITHARTIFACTLBL $WORK/unwarp/${id}_unwarp_phgseg_multilabel.nii.gz \
 -rm $AXISALIGN_SINGLE $WORK/unwarp/${id}_unwarp_phgsegshape_singlelabel.nii.gz \
 -rm $AXISALIGN_MULTI $WORK/unwarp/${id}_unwarp_phgsegshape_multilabel.nii.gz \
 -rm $AXISALIGN_SRLM $WORK/unwarp/${id}_unwarp_srlmseg_sr.nii.gz \
 -rm $AXISALIGN_SRLM_SHAPE $WORK/unwarp/${id}_unwarp_srlmsegshape_sr.nii.gz \
 -r $WORK/unwarp/${id}_unwarp.mat $CHAIN_SEG_TO_RAW 


 #Separate the multilabel image into separate components - csmedial and cslateral

  MULTILABEL_DIR=$WORK/unwarp/multi_label

  mkdir -p $MULTILABEL_DIR

  $C3D_HOME/c3d $WORK/unwarp/${id}_unwarp_phgseg_multilabel.nii.gz -retain-labels 3 4 \
 -replace 3 1 4 2 -o $MULTILABEL_DIR/${id}_unwarp_cslatseg.nii.gz

  $C3D_HOME/c3d $WORK/unwarp/${id}_unwarp_phgseg_multilabel.nii.gz -retain-labels 1 2 \
 -o $MULTILABEL_DIR/${id}_unwarp_csmedseg.nii.gz

  $C3D_HOME/c3d $WORK/unwarp/${id}_unwarp_phgseg_multilabel.nii.gz  -thresh 3 inf 1 0 \
 -o $MULTILABEL_DIR/${id}_unwarp_cslatsegshape.nii.gz

  $C3D_HOME/c3d $WORK/unwarp/${id}_unwarp_phgseg_multilabel.nii.gz -thresh 1 2 1 0 \
 -o $MULTILABEL_DIR/${id}_unwarp_csmedsegshape.nii.gz

 #Create separate label for the hippocampus, to help initialize registrations
 c3d -verbose $WORK/unwarp/${id}_unwarp_phgsegshape_singlelabel.nii.gz -as MTLSEG -stretch 0 1 -1 1 \
 -push MTLSEG $WORK/unwarp/${id}_unwarp_srlmseg_sr.nii.gz -int 0 -reslice-identity -thresh 1 1 -1 1 \
 -levelset 300 -thresh -inf 0 1 0 -push MTLSEG -add \
 -as X -thresh 1 1 1 0 -comp -thresh 1 1 1 0 -push X -thresh 2 2 1 0 \
 -o $WORK/unwarp/${id}_unwarp_hippo_roi_derived.nii.gz

 #Slightly widen the sulcus in the segmentations. Need to do this before applying any mesh smooting
 # For this to work well, need to ensure that there are no topological holes in the PHG segmentation
 #fn_widen_sulcus ${id}


}


fn_widen_sulcus()
{

 id=${1?}

 #Specify directories
 WDIR=$ROOT/preproc/${id}
 OUTPUT_DIR=$WDIR/sulcus_widen
 INPUT=$ROOT/preproc/${id}/${id}_axisalign_phgsegshape_singlelabel.nii.gz
 SRLM_SEG=$ROOT/preproc/${id}/${id}_axisalign_srlmseg_sr.nii.gz

 mkdir -p $OUTPUT_DIR

 #if [[ ! -f $OUTPUT_DIR/${id}_axisalign_srlm_mesh.stl ]]; then
 
  echo "Performing sulcus widening"

 # Get a binary mesh of the ROI
 $C3D_HOME/c3d $INPUT -trim 4vox -o $OUTPUT_DIR/binary_mesh.nii.gz

 # Compute hessian maps of sulcusness and gurusness
 $C3D_HOME/c3d $OUTPUT_DIR/binary_mesh.nii.gz -scale 1 -hessobj 2 0.2mm 0.3mm -o $OUTPUT_DIR/hessobj_gyr.nii.gz

 # Reduce the min-scale to 0.1mm
 $C3D_HOME/c3d $OUTPUT_DIR/binary_mesh.nii.gz -stretch 0 1 1 0 -hessobj 2 0.1mm 0.6mm -o $OUTPUT_DIR/hessobj_sulc.nii.gz

 # Compute erosion (opening sulci) and dilation (thickening cortex)
 $C3D_HOME/c3d $OUTPUT_DIR/binary_mesh.nii.gz -dilate 0 2x2x2 -o $OUTPUT_DIR/opening.nii.gz
 $C3D_HOME/c3d $OUTPUT_DIR/binary_mesh.nii.gz -dilate 1 2x2x2 -o $OUTPUT_DIR/closing.nii.gz

 # Perform voting between these two based on the hessian map
 $C3D_HOME/c3d $OUTPUT_DIR/hessobj_sulc.nii.gz -shift 0.000001 $OUTPUT_DIR/opening.nii.gz -stretch 0 1 1 0 -times \
    $OUTPUT_DIR/hessobj_gyr.nii.gz  -shift 0.000001 $OUTPUT_DIR/closing.nii.gz -times \
    -vote -comp -thresh 1 1 1 0 -o $OUTPUT_DIR/${id}_adjusted_mesh.nii.gz

 # Extract levelset
 vtklevelset $OUTPUT_DIR/${id}_adjusted_mesh.nii.gz $OUTPUT_DIR/${id}_axisalign_adjusted_mesh.stl 0.5

 #Extract levelset for srlm segmentation (don't need to adjust)
 vtklevelset $SRLM_SEG $OUTPUT_DIR/${id}_axisalign_srlm_mesh.stl 0.5

 #fi

}

## For atlas version 2, performed hippocampus label segmentation and cleaned up the MTL segmentation in a few places. Updating MTL segmentation and including green label for artifacts in hippocampus (saved in manual "inputs" directory")
function update_segmentations()
{
  # This is the subject ID
  id=${1?}

  ##Hippocampus label segmentation - contains medial/lateral extra hippocampal seg, hf seg, yellow and green artifacts.
 #Some green artifacts were overwritten by hf label so need to update this
  HF_SEG=$ROOT/inputs/$id/${id}_axisalign_hfseg_multilabel.nii.gz
  PHG_SEG=$ROOT/inputs/${id}/${id}_axisalign_phg_cs.nii.gz

  c3d $PHG_SEG -retain-labels 2 -as MED_ARTIFACT \
  $HF_SEG -thresh 1 6 1 0 -push MED_ARTIFACT -times -as MED_ARTIFACT_EDIT \
  $HF_SEG -replace 2 6 -push MED_ARTIFACT_EDIT -add -replace 5 1 7 2 8 2 -o $ROOT/inputs/${id}/${id}_axisalign_phg_cs_updated.nii.gz

}

function srlm_axisalign(){

  id=${1?}

  RAW_TO_DANSEG_MAT=$ATLAS_2016/preproc/${id}/${id}_raw_to_danseg.mat
  DANSEG_DB=$ATLAS_2016/input/${id}/${id}_dbseg.nii.gz
  RAW_TO_AXISALIGN_MAT=$ROOT/preproc/$id/${id}_transform_to_axisalign.mat

  AXISALIGN_IMG=$ROOT/inputs/$id/${id}_axisalign_img.nii.gz
  RAW_IMG=$ROOT/preproc/$id/${id}_raw_n4clip.nii.gz

 $GREEDY_HOME/greedy -d 3 \
 -rf $AXISALIGN_IMG \
 -ri LABEL 0.24vox \
 -rm $DANSEG_DB $ROOT/inputs/${id}/${id}_axisalign_danseg_db.nii.gz \
 -r $RAW_TO_AXISALIGN_MAT $ROOT/preproc/$id/unwarp/${id}_unwarp.mat $RAW_TO_DANSEG_MAT,-1

}


function main()
{
  mkdir -p $WDIR/dump


  N=$(cat $SUBJ_LIST | wc -l)

  for ((i=10;i<=10;i++)); do

 	id=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
	$ROOT/scripts/job_launcher.sh -m 25G -o $WDIR/dump -N "prep_atlas2022_${id}" $ROOT/scripts/organize_inputs.sh organize_subj_2022 $id

  # update_segmentations $id


  done


}

if [[ $1 ]]; then

  command=$1
  shift
  $command $@

else

  main
fi
