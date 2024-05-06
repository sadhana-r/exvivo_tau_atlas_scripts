#!/bin/bash
#$ -S /bin/bash
set -x -e

# The directory where we are doing this work
ROOT=/project/hippogang_3/sravikumar/atlasPHG2019

# The work directory
WDIR=$ROOT/atlasPHG_V3_2022/mst_multilabel

# Read the common include file
source $ROOT/scripts/common.sh

# List of all the ids to run
SUBJ_LIST=$ROOT/scripts/subj_atlas_2022_v2.txt
IDS=$(cat $SUBJ_LIST)

# Global parameters - reused by script - initially all 5
GSHOOT_NITER=5
GSHOOT_NITER_NCC=5
NCC_NITER=6
NSLOTS=4

TMPDIR=/tmp

# Function : return the input image, phgseg for a subject
# use with read a b <<<$(fn_input_unwarped_images $id)
function fn_input_unwarped_images()
{
  local ID=${1?}

  local SEG_SHAPEONLY=${2?}

  # Template is being built from images that have been unwarped to correct
  # for scanner gradient scaling

  if [[ $SEG_SHAPEONLY -eq 1 ]]; then
	  local IM=$ROOT/preproc/$ID/unwarp/${ID}_unwarp_img.nii.gz
	  local CSMED=$ROOT/preproc/$ID/unwarp/multi_label/${ID}_unwarp_csmedsegshape.nii.gz
          local CSLAT=$ROOT/preproc/$ID/unwarp/multi_label/${ID}_unwarp_cslatsegshape.nii.gz
	  local PHG=$ROOT/preproc/$ID/unwarp/${ID}_unwarp_phgsegshape_singlelabel.nii.gz
 	  local SRLM=$ROOT/preproc/$ID/unwarp/${ID}_unwarp_srlmsegshape_sr.nii.gz

	  echo $IM $CSMED $CSLAT $PHG $SRLM

# This segmentation has a separate label marking for regions where the intensity information should be ignored
else

          local IM=$ROOT/preproc/$ID/unwarp/${ID}_unwarp_img.nii.gz
          local CSMED=$ROOT/preproc/$ID/unwarp/multi_label/${ID}_unwarp_csmedseg.nii.gz
          local CSLAT=$ROOT/preproc/$ID/unwarp/multi_label/${ID}_unwarp_cslatseg.nii.gz
	  local PHG=$ROOT/preproc/$ID/unwarp/${ID}_unwarp_phgseg_singlelabel.nii.gz
          local SRLM=$ROOT/preproc/$ID/unwarp/${ID}_unwarp_srlmseg_sr.nii.gz

	  echo $IM $CSMED $CSLAT $PHG $SRLM

fi
}


function fn_all_pairwise()
{
  # Launch all jobs
  N=$(cat $SUBJ_LIST | wc -l)

  for ((i=1;i<=$N;i++)); do
    ID1=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
    for ((j=31;j<=$N;j++)); do
      if [[ $i -ne $j ]]; then
        ID2=$(cat $SUBJ_LIST | head -n $j | tail -n 1)
           $ROOT/scripts/job_launcher.sh -m 15G -o $WDIR/dump -N "pw_${ID1}_${ID2}" \
          "$0" fn_pairwise_rough $ID1 $ID2
      fi
    done
  done

  # Wait for completion
  $ROOT/scripts/job_launcher.sh -w "pw_*" /bin/sleep 0
}

# Function : perform 'smooth' registration between pair of scans in order
# to evaluate similarity
function fn_pairwise_rough()
{
  # The ids of the two subjects
  IDREF=${1?}
  IDMOV=${2?}

  SEG_SHAPEONLY=1;

  # The images being used
  read IMGREF CSMEDREF CSLATREF PHGREF SRLMREF <<<$(fn_input_unwarped_images $IDREF $SEG_SHAPEONLY)
  read IMGMOV CSMEDMOV CSLATMOV PHGMOV SRLMMOV <<<$(fn_input_unwarped_images $IDMOV $SEG_SHAPEONLY)

  #Input hippomcapus ROI segmentations
  HIPPO_REF=$ROOT/preproc/$IDREF/unwarp/${IDREF}_unwarp_hippo_roi_derived.nii.gz
  HIPPO_MOV=$ROOT/preproc/$IDMOV/unwarp/${IDMOV}_unwarp_hippo_roi_derived.nii.gz

  # The tag for this pair
  TAG="ref_${IDREF}_mov_${IDMOV}"

  # Registration will be performed on the label masks of the two images
  # using a very high smoothness term
  WORK=$WDIR/pairwise/${TAG}
  mkdir -p $WORK

  # warp file
  MAT_MOMENTS=$WORK/${TAG}_moments.mat
  MAT_RIGID=$WORK/${TAG}_rigid.mat
  MAT_AFFINE=$WORK/${TAG}_affine.mat
  WARP=$WORK/${TAG}_warp.nii.gz

  # Resliced moving segmentations
  RSLPATTERN=$WORK/${TAG}_reslice_XXXseg.nii.gz
  CSMEDRSL=$WORK/${TAG}_reslice_csmedseg.nii.gz
  CSLATRSL=$WORK/${TAG}_reslice_cslatseg.nii.gz
  SRLMRSL=$WORK/${TAG}_reslice_srlmseg.nii.gz


  for fn in $CSMEDREF $CSLATREF $SRLMREF $CSMEDMOV $CSLATMOV $SRLMMOV; do
    ln -sf $fn $WORK/
  done

  # This line of greedy inputs and weights will be reused - was 1 1 5
  GREEDY_INPUTS="-w 1 -i $CSMEDREF $CSMEDMOV -w 2 -i $CSLATREF $CSLATMOV -w 5 -i $SRLMREF $SRLMMOV"

<<"IGNORE"
#Works better wth hippocampus segmentations. Don't use sides - not reliable
  # Do we need to flip space between these two segmentations?
  SIDE_REF=$(echo $IDREF | awk -F '-' '{print $NF}')
  SIDE_MOV=$(echo $IDMOV | awk -F '-' '{print $NF}')

  if [[ $SIDE_REF == $SIDE_MOV ]]; then DET="1"; else DET="-1"; fi
IGNORE

if [[ ! -f $WARP ]] ; then

 # Perform moments of intertia matching between the two masks - added NMI metric
 # Moment matching outputs NaN when I include both inputs. Using only medial segmentation
 # to perform moment matching. Affine will use both.
 # Aug 2022- moments on CSMED doesn;t seem to work. MTL often flipped. Try SRLM?
 $GREEDY_HOME/greedy -d 3 -threads $NSLOTS \
 -i $HIPPO_REF $HIPPO_MOV \
 -moments \
 -o $MAT_MOMENTS

<<"NOTWORKING"
#remove -det $DET
  $GREEDY_HOME/greedy -d 3 -a -threads $NSLOTS \
 -i $HIPPO_REF $HIPPO_MOV \
 -ia-image-centers \
 -dof 6 -o $MAT_RIGID -n 100x100x0 -m NCC 4x4x4 -search 2000 any 20
NOTWORKING

 # Perform affine matching between the two masks
 $GREEDY_HOME/greedy -d 3 -a \
 $GREEDY_INPUTS \
 -ia $MAT_MOMENTS \
 -n 100x100 \
 -o $MAT_AFFINE

  #Reduce iterations from 100/80/80 to 100/80/40
  # Run greedy between these two images - included NCC metric
  $GREEDY_HOME/greedy -d 3 -threads $NSLOTS \
  $GREEDY_INPUTS \
  -it $MAT_AFFINE \
  -n 100x80x40x0 \
  -s 2.0mm 0.1mm -e 0.5 -m NCC 4x4x4 \
  -o $WARP

  # Reslice the segmentations from raw space
  reslice_subj_segshapesonly $IDMOV $CSMEDREF $RSLPATTERN $WARP $MAT_AFFINE

  # Compute the overlaps - 1 = csmed, 2 = cslat, 3 = srlm
  $C3D_HOME/c3d \
  $CSLATREF -scale 2 -as LAT $SRLMREF -scale 2 $CSMEDREF -add -push LAT -add -popas REF \
  $CSLATRSL -scale 2 -as LAT $SRLMRSL -scale 2 $CSMEDRSL -add -push LAT -add \
  -push REF -label-overlap | tee $WORK/${TAG}_overlap.txt

  # Write just the adjacency value
  cat $WORK/${TAG}_overlap.txt | awk 'NR==3 {print $3}' > $WORK/${TAG}_adj.txt
  cat $WORK/${TAG}_overlap.txt | awk 'NR==6 {print $4}' > $WORK/${TAG}_adj_1.txt
  cat $WORK/${TAG}_overlap.txt | awk 'NR==7 {print $4}' > $WORK/${TAG}_adj_2.txt
  cat $WORK/${TAG}_overlap.txt | awk 'NR==8 {print $4}' > $WORK/${TAG}_adj_3.txt
fi

}

# Private function used by reslice_subj() and reslice_subj_segsonly()
function reslice_subj_internal()
{
  # Mode flags
  local DO_IMG=${1?}
  local DO_SEG=${2?}
  local SEG_SHAPEONLY=${3?}

# Subject's ID
  local ID=${4?}

  # Reference space
  local REFSPACE=${5?}
  local OUTPATTERN=${6?}

  # The requested warp chain
  shift 6
  local WARPCHAIN="$@"

  # Subject's preproc directory
  local PDIR=$ROOT/preproc/${ID}

  # Subject's RAW MRI scan
  local RAW_IMG=$PDIR/${ID}_raw_n4clip.nii.gz

  # Subject's source segmentation of PHG
  if [[ $SEG_SHAPEONLY -eq 1 ]]; then
  	local CSMEDSEG=$PDIR/${ID}_axisalign_csmedsegshape.nii.gz
	local CSLATSEG=$PDIR/${ID}_axisalign_cslatsegshape.nii.gz
	local PHGSEG=$PDIR/${ID}_axisalign_phgsegshape_singlelabel.nii.gz
	local SRLMSEG=$PDIR/${ID}_axisalign_srlmseg_sr.nii.gz
  else
	local CSMEDSEG=$PDIR/${ID}_axisalign_csmedseg.nii.gz
	local CSLATSEG=$PDIR/${ID}_axisalign_cslatseg.nii.gz
	local PHGSEG=$PDIR/${ID}_axisalign_phgseg_singlelabel.nii.gz
	local SRLMSEG=$PDIR/${ID}_axisalign_srlmseg_sr.nii.gz
  fi

  # Subject's unwarp and axisalign/raw matrices
  local UNWARP_MAT=$PDIR/unwarp/${ID}_unwarp.mat

<<"DONTNEED"
  # Extract format of axisalign to raw warpchain - saved in seg_to_raw.txt
  SEG_RAW_FORMAT=$(cat $ROOT/scripts/seg_to_raw.txt | awk -v id="$ID" '$1 == id {print $2}')

 if [[ $SEG_RAW_FORMAT -eq 1 ]]; then
	local SEG_TO_RAW_MAT="$PDIR/${ID}_warp.mat $PDIR/${ID}_inverse_axisalign.mat"
 fi

 if [[ $SEG_RAW_FORMAT -eq 2 ]]; then
	local SEG_TO_RAW_MAT="$PDIR/${ID}_inverse_axisalign.mat"
 fi
DONTNEED

  SEG_TO_RAW_MAT="$PDIR/${ID}_transform_to_axisalign.mat,-1"

  # Warp chain from raw image space to template space
  local CHAIN_RAW_TO_TEMPLATE="$WARPCHAIN $UNWARP_MAT"

  # Warp chain from input image space (where segmentation drawn) to template
  local CHAIN_INPUT_TO_TEMPLATE="$CHAIN_RAW_TO_TEMPLATE $SEG_TO_RAW_MAT"

  # Reslice the image
  if [[ $DO_IMG -eq 1 ]]; then
    $GREEDY_HOME/greedy -d 3 \
      -rf $REFSPACE \
      -rm $RAW_IMG $(echo $OUTPATTERN | sed -e "s/XXX/img/g") \
      -r $CHAIN_RAW_TO_TEMPLATE
  fi

  # Reslice the segmentations
  if [[ $DO_SEG -eq 1 ]]; then

	if [[ $SEG_SHAPEONLY -eq 1 ]]; then

    		$GREEDY_HOME/greedy -d 3 \
      		-rf $REFSPACE \
      		-ri LABEL 0.24vox \
      		-rm $CSMEDSEG $(echo $OUTPATTERN | sed -e "s/XXX/csmed/g") \
                -rm $CSLATSEG $(echo $OUTPATTERN | sed -e "s/XXX/cslat/g") \
		-rm $PHGSEG $(echo $OUTPATTERN | sed -e "s/XXX/phg/g") \
		-rm $SRLMSEG $(echo $OUTPATTERN | sed -e "s/XXX/srlm/g") \
      		-r $CHAIN_INPUT_TO_TEMPLATE
	else
		$GREEDY_HOME/greedy -d 3 \
                -rf $REFSPACE \
                -ri LABEL 0.24vox \
                -rm $CSMEDSEG $(echo $OUTPATTERN | sed -e "s/XXX/csmed_int/g") \
		-rm $CSLATSEG $(echo $OUTPATTERN | sed -e "s/XXX/cslat_int/g") \
		-rm $PHGSEG $(echo $OUTPATTERN | sed -e "s/XXX/phg_int/g") \
		-rm $SRLMSEG $(echo $OUTPATTERN | sed -e "s/XXX/srlm/g") \
                -r $CHAIN_INPUT_TO_TEMPLATE

	fi
  fi
}

# This function reslices the subject's RAW image data and native space
# segmentations into the provided reference image space using the chain
# of transforms between the reference space and the subject's UNWARPED
# space.
#
# Usage:
#   reslice_subj ID reference output_pattern warp_chain
#
# Output pattern must include string XXX which will be replaced by
# img, phg  respectively.
function reslice_subj()
{
  reslice_subj_internal 1 1 1 $@
}


# Same as above but segmentations only
function reslice_subj_segsonly()
{
  reslice_subj_internal 0 1 0 $@
}

# Same as above but segmentation shapes only
function reslice_subj_segshapesonly()
{
  reslice_subj_internal 0 1 1 $@
}

function fn_print_adjacency()
{

  # Create the adjacency matrix
  N=$(cat $SUBJ_LIST | wc -l)

<<"SKIP"
  LABELS="1 2 3"
  for label in $LABELS; do
  	for ((i=1;i<=$N;i++)); do
    	IDi=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
    	  for ((j=1;j<=$N;j++)); do
    	  	if [[ $i -ne $j ]]; then
        	  IDj=$(cat $SUBJ_LIST | head -n $j | tail -n 1)
        	  TAG="ref_${IDi}_mov_${IDj}"
        	  if [[ -f $WDIR/pairwise/$TAG/${TAG}_adj_${label}.txt ]]; then
        		  echo $i $j $(cat $WDIR/pairwise/$TAG/${TAG}_adj_${label}.txt)
        	  fi
      		fi
    	  done
  	done > $WDIR/adj_${label}.txt
   done
SKIP
   
for ((i=1;i<=$N;i++)); do
        IDi=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
          for ((j=1;j<=$N;j++)); do
                if [[ $i -ne $j ]]; then
                  IDj=$(cat $SUBJ_LIST | head -n $j | tail -n 1)
                  TAG="ref_${IDi}_mov_${IDj}"
                  if [[ -f $WDIR/pairwise/$TAG/${TAG}_adj.txt ]]; then
                          echo $i $j $(cat $WDIR/pairwise/$TAG/${TAG}_adj.txt)
                  fi
                fi
          done
   done > $WDIR/adj.txt

# Use R to compute paths
Rscript $ROOT/scripts/compute_mst.R $WDIR/adj.txt > $WDIR/paths.txt

# Rscript $ROOT/scripts/compute_mst.R $WDIR/adj_1.txt  $WDIR/adj_2.txt $WDIR/adj_3.txt > $WDIR/paths.txt


}

function fn_register_root_to_oldatlas()
{
  # We would like the MST root to be correctly aligned with the old atlas so
  # it's not strangely rotated.
  local ROOTIDX=$(cat $WDIR/paths.txt | awk 'NF==1 {print $1}')
  local ROOTID=$(cat $SUBJ_LIST | head -n $ROOTIDX | tail -n 1)

  #Find the transformation from the raw image space to the axisaligned space (instead of using the piecewise transformations)
  SEG_SHAPEONLY=1
  AXISALIGN_CSMED=$ROOT/preproc/${ROOTID}/${ROOTID}_axisalign_csmedsegshape.nii.gz
  AXISALIGN_CSLAT=$ROOT/preproc/${ROOTID}/${ROOTID}_axisalign_cslatsegshape.nii.gz
  AXISALIGN_SRLM=$ROOT/preproc/${ROOTID}/${ROOTID}_axisalign_srlmseg_sr.nii.gz

 # The output directory
  local WORK=$WDIR/root_to_oldatlas

  mkdir -p $WORK

  #Create a symbolic link to the axisaligned segmentation for use by other functions
  ln -sf $AXISALIGN_CSMED $WORK/root_axisalign_csmedsegshape.nii.gz
  ln -sf $AXISALIGN_CSLAT $WORK/root_axisalign_cslatsegshape.nii.gz
  ln -sf $AXISALIGN_SRLM $WORK/root_axisalign_srlmseg.nii.gz

  # The input hippo segmentation
  local ROOT_IMG ROOT_CSMED ROOT_CSLAT ROOT_PHG ROOT_SRLM
  read ROOT_IMG ROOT_CSMED ROOT_CSLAT ROOT_PHG ROOT_SRLM  <<<$(fn_input_unwarped_images $ROOTID $SEG_SHAPEONLY)

  # Added this because when the root id 120267, the unwarped segmentation is right by the image border.
  # Registrations are getting cut off.
  ROOT_IMG_PAD=$WORK/root_img_pad.nii.gz
  ROOT_CSMED_PAD=$WORK/root_csmed_pad.nii.gz
  ROOT_CSLAT_PAD=$WORK/root_cslat_pad.nii.gz
  ROOT_SRLM_PAD=$WORK/root_srlm_pad.nii.gz

<<"120267"
  $C3D_HOME/c3d $ROOT_IMG -pad 30x0x0 0x0x0 0 -o $ROOT_IMG_PAD
  $C3D_HOME/c3d $ROOT_CSMED -pad 30x0x0 0x0x0 0 -o $ROOT_CSMED_PAD
  $C3D_HOME/c3d $ROOT_CSLAT -pad 30x0x0 0x0x0 0 -o $ROOT_CSLAT_PAD
  $C3D_HOME/c3d $ROOT_SRLM -pad 30x0x0 0x0x0 0 -o $ROOT_SRLM_PAD
120267

  #When root is HNL39 - need to padd in the other direction. So instead, just pad 10 in each dir.
  $C3D_HOME/c3d $ROOT_IMG -pad 10x10x10 30x10x10 0 -o $ROOT_IMG_PAD
  $C3D_HOME/c3d $ROOT_CSMED -pad 10x10x10 30x10x10 0 -o $ROOT_CSMED_PAD
  $C3D_HOME/c3d $ROOT_CSLAT -pad 10x10x10 30x10x10 0 -o $ROOT_CSLAT_PAD
  $C3D_HOME/c3d $ROOT_SRLM -pad 10x10x10 30x10x10 0 -o $ROOT_SRLM_PAD


  # The matrices output
  local MAT_MOMENTS=$WORK/root_to_atlas_moments.mat
  local MAT_RIGID=$WORK/root_to_atlas_rigid.mat

  # Do the registration
  # Perform moments of intertia matching between the two masks
  $GREEDY_HOME/greedy -d 3 \
    -i $AXISALIGN_CSMED $ROOT_CSMED_PAD \
    -i $AXISALIGN_CSLAT $ROOT_CSLAT_PAD \
    -moments -o $MAT_MOMENTS

  # Perform affine matching between the two masks
  $GREEDY_HOME/greedy -d 3 -threads $NSLOTS \
    -i $AXISALIGN_CSMED $ROOT_CSMED_PAD \
    -i $AXISALIGN_CSLAT $ROOT_CSLAT_PAD \
    -a -dof 6 -ia $MAT_MOMENTS \
    -n 100x100 \
    -o $MAT_RIGID
}

function fn_all_paths()
{
  
 func=${1?}
 if [[ -f $WDIR/paths ]]; then
    rm -r $WDIR/paths
  fi
 
  # Launch all jobs
  N=$(cat $SUBJ_LIST | wc -l)
  for ((i=1;i<=$N;i++)); do
 	$ROOT/scripts/job_launcher.sh -m 15G -o $WDIR/dump -N "path_$i" \
          "$0" $func $i
  done


  # Wait for completion
  $ROOT/scripts/job_launcher.sh -w "path_*" /bin/sleep 0
}

function fn_register_mst()
{
  # The index of the subject being registered using the MST
  idx=${1?}

  # The ID of the subject
  ID=$(cat $SUBJ_LIST | head -n $idx | tail -n 1)

  SEG_SHAPEONLY=1

  # The path of nodes leading up to the root node (from current to the root)
  GPATH=($(cat $WDIR/paths.txt | head -n $idx | tail -n 1))

  # The native space images and registrations
  read IMGNAT CSMEDNAT CSLATNAT PHGNAT SRLMNAT <<<$(fn_input_unwarped_images $ID $SEG_SHAPEONLY)

  # We start with the warp chain as an empty list
  WARPCHAIN=""

  # The moving images at the start are the same as the native images
  IMGMOV=$IMGNAT
  CSMEDMOV=$CSMEDNAT
  CSLATMOV=$CSLATNAT
  SRLMMOV=$SRLMNAT

  # The id in whose space the moving image is - at the beginning this is the
  # same as the subject's ID
  IDMOV=$ID

  # Loop. The first element of GPATH is the image itself
  for ((j=1;j<${#GPATH[*]};++j)); do

    # The next image in the path
    idxref=${GPATH[j]}
    IDREF=$(cat $SUBJ_LIST | head -n $idxref | tail -n 1)

    # The target images
    read IMGREF CSMEDREF CSLATREF PHGREF SRLMREF <<<$(fn_input_unwarped_images $IDREF $SEG_SHAPEONLY)

    # Create directory for this registration
    WORK=$WDIR/paths/${ID}/step_${j}
    mkdir -p $WORK

    # Added this because when the root id 120267, the unwarped segmentation is right by the image border.
    # Registrations are getting cut off.
    IMGREF_PAD=$WORK/img_ref_pad.nii.gz
    CSMEDREF_PAD=$WORK/csmed_ref_pad.nii.gz
    CSLATREF_PAD=$WORK/cslat_ref_pad.nii.gz
    SRLMREF_PAD=$WORK/srlm_ref_pad.nii.gz

  # when root is HNL39, need to pad in other direction. Added 10 vox to be safe on other sides
    $C3D_HOME/c3d $IMGREF -pad 10x10x10 30x10x10 0 -o $IMGREF_PAD
    $C3D_HOME/c3d $CSMEDREF -pad 10x10x10 30x10x10 0 -o $CSMEDREF_PAD
    $C3D_HOME/c3d $CSLATREF -pad 10x10x10 30x10x10 0 -o $CSLATREF_PAD
    $C3D_HOME/c3d $SRLMREF -pad 10x10x10 30x10x10 0 -o $SRLMREF_PAD

    # We already know the affine between these two images
    AFFINE_INIT=$WDIR/pairwise/ref_${IDREF}_mov_${IDMOV}/ref_${IDREF}_mov_${IDMOV}_affine.mat

    # The output warp
    AFFINE=$WORK/affine.mat
    WARP=$WORK/warp.nii.gz

    # Create links for the relevant images
    ln -sf $IMGMOV $WORK/img_moving.nii.gz
    ln -sf $CSMEDMOV $WORK/csmed_moving.nii.gz
    ln -sf $CSLATMOV $WORK/cslat_moving.nii.gz
    ln -sf $SRLMMOV $WORK/srlm_moving.nii.gz

    ln -sf $IMGREF $WORK/img_ref.nii.gz
    ln -sf $CSMEDREF $WORK/csmed_ref.nii.gz
    ln -sf $CSLATREF $WORK/cslat_ref.nii.gz
    ln -sf $SRLMREF $WORK/srlm_ref.nii.gz

    ln -sf $AFFINE_INIT $WORK/affine_init.mat

  # Perform affine registration between the moving and fixed
    $GREEDY_HOME/greedy -d 3 \
      -w 1000 -i $CSMEDREF_PAD $CSMEDMOV \
      -w 2000 -i $CSLATREF_PAD $CSLATMOV \
      -w 5000 -i $SRLMREF_PAD $SRLMMOV \
      -threads $NSLOTS -n 100x100 \
      -ia $AFFINE_INIT -a -o $AFFINE

    # Perform registration - use intensity along with boundary strength - took out -w 0.00006 -i $IMGREF $IMGMOV
	# Increase number of iterations from 100x50x40x20. Chnage interations from 100/80/50/30 to 100/80/40/20
      $GREEDY_HOME/greedy -d 3 \
      -w 1000 -i $CSMEDREF_PAD $CSMEDMOV \
      -w 2000 -i $CSLATREF_PAD $CSLATMOV \
      -w 5000 -i $SRLMREF_PAD $SRLMMOV \
      -w 0.0007 -i $IMGREF_PAD $IMGMOV \
      -threads $NSLOTS -n 100x80x40x20 \
      -s 0.6mm 0.1mm -e 0.5 -float \
      -it $AFFINE \
      -o $WARP

    # Update the warp chain and the reslices
    WARPCHAIN="$WARP $AFFINE $WARPCHAIN"
    RESLICE_PTRN=$WORK/reslice_${ID}_to_${IDREF}_XXX.nii.gz
    IMGMOV=$WORK/reslice_${ID}_to_${IDREF}_img.nii.gz
    CSMEDMOV=$WORK/reslice_${ID}_to_${IDREF}_csmed.nii.gz
    CSLATMOV=$WORK/reslice_${ID}_to_${IDREF}_cslat.nii.gz
    SRLMMOV=$WORK/reslice_${ID}_to_${IDREF}_srlm.nii.gz
    IDMOV=${IDREF}

    # Apply the deformations to the data - going all the way to the raw IMG and SEG - this is using
    # the segmenation shapes only
    reslice_subj $ID $CSMEDREF_PAD $RESLICE_PTRN $WARPCHAIN

 done

  # In the final stage, we prepend the warpchain by the rigid trasform to the
  # old atlas reference space, so that the atlases live in similar space and
  # so that the orientation of the hippocampus is reasonable
  WORK=$WDIR/paths/${ID}/final

  # The final reference space comes from the old hippocampus atlas
  WARPCHAIN="$WDIR/root_to_oldatlas/root_to_atlas_rigid.mat $WARPCHAIN"
  RESLICE_PTRN=$WORK/reslice_${ID}_to_${IDREF}_XXX.nii.gz
  AXISALIGN_SPACE=$WDIR/root_to_oldatlas/root_axisalign_csmedsegshape.nii.gz

  # Create links to the "final" moving images
  mkdir -p $WORK
  reslice_subj $ID $AXISALIGN_SPACE $RESLICE_PTRN $WARPCHAIN

  # Save the warpchain to a text file for easy reuse
  echo $WARPCHAIN > $WORK/chain_unwarp_to_final.txt

}

function fn_paths_average()
{
  # Create the mean image and mean DB images
  mkdir -p $WDIR/template_mst
  for what in img csmed cslat srlm phg; do
    $C3D_HOME/c3d $WDIR/paths/*/final/reslice_*${what}.nii.gz \
      -mean -o $WDIR/template_mst/template_mst_${what}.nii.gz
  done
}

function fn_pair_paths_overlap()
{
  ID1=${1?}
  ID2=${2?}

  WORK=$WDIR/template_mst/overlap/pairs
  mkdir -p $WORK

  # Compute Dice
  for what in phg srlm; do
    $C3D_HOME/c3d \
      $WDIR/paths/${ID1}/final/*_${what}*.nii.gz \
      $WDIR/paths/${ID2}/final/*_${what}*.nii.gz \
      -overlap 1 > $WORK/overlap_${what}_${ID1}_${ID2}.txt
  done

  # Compute NCC
  c3d \
    $WDIR/paths/${ID1}/final/*_img*.nii.gz \
    $WDIR/paths/${ID2}/final/*_img*.nii.gz \
    -ncc 4x4x4 -replace -NaN 0 NaN 0 -Inf 0 Inf 0 \
    $WDIR/template_mst/template_mst_phg.nii.gz -thresh 0.5 inf 1 0 \
    -lstat | awk 'NR==3 {print $2}' \
    > $WORK/ncc_${ID1}_${ID2}.txt

}

function fn_all_pair_paths_overlap()
{
  func=${1?}

  N=$(cat $SUBJ_LIST | wc -l)
  for ((i=1;i<=$N;i++)); do
    IDi=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
    for ((j=$((i+1));j<=$N;j++)); do
      if [[ $i -ne $j ]]; then
        IDj=$(cat $SUBJ_LIST | head -n $j | tail -n 1)
        $ROOT/scripts/job_launcher.sh -m 150G -o $WDIR/dump -N "ovl_${IDi}_${IDj}" \
          "$0" $func $IDi $IDj
      fi
    done
  done

  # Wait for completion
  $ROOT/scripts/job_launcher.sh -w "ovl_*" /bin/sleep 0
}

function fn_shooting_init_landmarks()
{
  # Determine what is the input
  MODE=${1?}
  if [[ $MODE == 'MST' ]]; then
    IMG_CSMED=$WDIR/template_mst/template_mst_csmed.nii.gz
    IMG_CSLAT=$WDIR/template_mst/template_mst_cslat.nii.gz
    IMG_SRLM=$WDIR/template_mst/template_mst_srlm.nii.gz
    WORK=$WDIR/gshoot/template
  elif [[ $MODE == 'NCC' ]]; then
    IMG_CSMED=$(ls $WDIR/template_ncc/iter*/template_csmed.nii.gz | tail -n 1)
    IMG_CSLAT=$(ls $WDIR/template_ncc/iter*/template_cslat.nii.gz | tail -n 1)
    IMG_SRLM=$(ls $WDIR/template_ncc/iter*/template_srlm.nii.gz | tail -n 1)
    WORK=$WDIR/gshoot_ncc/template
  fi

  # Samples
  SAMCSMED=$WORK/work/samples_csmed.ply
  SAMCSLAT=$WORK/work/samples_cslat.ply
  SAMSRLM=$WORK/work/samples_srlm.ply
  TEMPLATE=$WORK/root_landmarks.vtk

  # Create output directories
  mkdir -p $WORK/work

  # Extract the level set surfaces
  vtklevelset $IMG_CSMED $WORK/work/mesh_csmed.stl 0.5
  vtklevelset $IMG_CSLAT $WORK/work/mesh_cslat.stl 0.5
  vtklevelset $IMG_SRLM $WORK/work/mesh_srlm.stl 0.5

  # Perform point sampling
  /project/hippogang_2/pauly/bin/mesh_poisson_sample $WORK/work/mesh_csmed.stl $SAMCSMED 500
  /project/hippogang_2/pauly/bin/mesh_poisson_sample $WORK/work/mesh_cslat.stl $SAMCSLAT 250
  /project/hippogang_2/pauly/bin/mesh_poisson_sample $WORK/work/mesh_srlm.stl $SAMSRLM 500

  # Combine the two PLY objects into a VTK mesh
  NVCSMED=$(cat $SAMCSMED | grep 'element vertex' | awk '{print $3}')
  NVCSLAT=$(cat $SAMCSLAT | grep 'element vertex' | awk '{print $3}')
  NVSRLM=$(cat $SAMSRLM | grep 'element vertex' | awk '{print $3}')
  NV=$((NVCSMED + NVCSLAT + NVSRLM))

# Write the header of the VTK file
  echo "# vtk DataFile Version 4.0" > $TEMPLATE
  echo "vtk output" >> $TEMPLATE
  echo "ASCII" >> $TEMPLATE
  echo "DATASET POLYDATA" >> $TEMPLATE
  echo "POINTS ${NV} float" >> $TEMPLATE

  # Write the point coordinates
  grep -A $NVCSMED end_header $SAMCSMED | tail -n $NVCSMED >> $TEMPLATE
  grep -A $NVCSLAT end_header $SAMCSLAT | tail -n $NVCSLAT >> $TEMPLATE
  grep -A $NVSRLM end_header $SAMSRLM | tail -n $NVSRLM >> $TEMPLATE

}

function fn_shooting_correction_subj_iter()
{
  # The index of the subject being registered using the MST
  idx=${1?}

  # The iteration - 0 is the initial iteration, involves extra work
  iter=${2?}

  # The path to the landmarks
  LANDMARKS=${3?}

  # Which mode are we running shooting in?
  # MST - after building the initial MST
  # NCC - after building the image-based template
  MODE=${4?}

 # Use the segmentation shapes only at this stage
  SEG_SHAPEONLY=1

  # The ID of the subject
  ID=$(cat $SUBJ_LIST | head -n $idx | tail -n 1)

  # Output directory - depends on the mode
  if [[ $MODE == 'MST' ]]; then
    WORK=$WDIR/gshoot/${ID}
  elif [[ $MODE == 'NCC' ]]; then
    WORK=$WDIR/gshoot_ncc/${ID}
  fi

  WITER=$WORK/iter_${iter}
  mkdir -p $WITER $WITER/movie

  # The native space images and registrations
  read IMGNAT CSMEDNAT CSLATNAT PHGNAT SRLMNAT <<<$(fn_input_unwarped_images $ID $SEG_SHAPEONLY)

  # Reference space (all subjects resliced to the old atlas space (axis align))
  REFSPACE=$WORK/refspace.nii.gz

  # Result meshes
  TARGET=$WORK/shooting_target_native.vtk
  LM_PROCRUSTES=$WITER/shooting_target_procrustes.vtk
  LM_PROCRUSTES_MAT=$WITER/target_to_root_procrustes.mat
  SHOOTING_WARP=$WITER/shooting_warp.nii.gz

  # Output pattern
  RESLICE_PTRN=$WITER/reslice_${ID}_shooting_to_template_XXX.nii.gz

  # Mesh containing the momenta
  MOMENTA=$WITER/shooting_momenta.vtk

  #if [[ ! -f $SHOOTING_WARP ]]; then
  # Target-related stuff in the WORK directory that is only done in the first iter
  if [[ $iter -eq 0 ]]; then

    # Create refspace - the funny wildcard here makes sure that the reference
    # space is also picked up for the root subject in the MST
    REFSPACE_SRC=$(ls $WDIR/paths/${ID}/final/*_csmed*.nii.gz)
    $C3D_HOME/c3d $REFSPACE_SRC -pad 0x20x20 0x20x20 0 -o $REFSPACE

# Compute the target landmarks for this subject
    if [[ $MODE == 'MST' ]]; then

      # Get the warp chain from file
      WARPCHAIN=$(cat $WDIR/paths/${ID}/final/chain_unwarp_to_final.txt)

      # Apply the warp chain to the landmark mesh in template space, creating
      # the target locations for the geodesic shooting registration
      # -- this code works when the WARPCHAIN is empty (MST root).
      # Obtain the landmarks for each specimen in its native space by undoing the MST warp
      # TARGET is the landmarks in the native (unwarp) space and LANDMARKS is the landmarks in the root space
      $GREEDY_HOME/greedy -d 3 \
        -rf $REFSPACE \
        -rs $LANDMARKS $TARGET \
        -r $WARPCHAIN

    elif [[ $MODE == 'NCC' ]]; then

      # The output of the segmentation-guided shooting iteration
      GSDIR=$WDIR/gshoot/${ID}/iter_$((GSHOOT_NITER-1))

      # The final NCC warp
      NCCWARP="$WDIR/template_ncc/iter$((NCC_NITER-1))/warp_root_template_to_${ID}.nii.gz,-64"

      # Build up the warp chain - all the 'shape' warps
      WARPCHAIN="$NCCWARP $GSDIR/shooting_warp.nii.gz $GSDIR/target_to_root_procrustes.mat,-1"

      # Apply the warp chain to the landmark mesh in template space, creating
      # the target locations for the geodesic shooting registration
      $GREEDY_HOME/greedy -d 3 \
        -rf $REFSPACE \
        -rs $LANDMARKS $TARGET \
        -r $WARPCHAIN

    fi

  fi

  # Landmarks in reference space
  ln -sf $LANDMARKS $WITER/landmarks.vtk

  # Bring the target mesh back near the root mesh using procrustes alignment - computes LM_PROCRUSTES_MAT
  vtkprocrustes $TARGET $LANDMARKS $LM_PROCRUSTES_MAT

  # Apply procrustes to the landmarks. Warp mesh is more reliable than greedy around
  # image edges, so we use it here. We had a situation in 65240-R where some landmarks
  # that were close to the image border were mapped badly (to a far-off location) which
  # screwed up the geodesic (and might have impacted the rest of the template too...
  # Applies LM_PROCRUSTES_MAT to move the native space landmarks to bring them closer to the root?
  warpmesh $TARGET $LM_PROCRUSTES $LM_PROCRUSTES_MAT

  # Perform geodesic shooting between the procrustes landmarks and the
  # warped landmarks - this is going to allow us to interpolate the correspondence
  # found by the MST to the rest of the images
  
time lmshoot -d 3 \
  -m $LANDMARKS $LM_PROCRUSTES \
  -o $MOMENTA \
  -s 2.0 -l 5000 -n 40 -i 240 0 -f -O $WITER/movie/movie%04d.vtk

# Convert the shooting result into a warp
  lmtowarp -d 3 -n 40 -r $REFSPACE \
   -m $MOMENTA -o $SHOOTING_WARP \
   -s 2.0


 if [[ $iter -eq $((GSHOOT_NITER-1)) ]]; then

    MOMENTA_INVERSE=$WITER/shooting_inverse_momenta.vtk
    SHOOTING_INVWARP=$WITER/shooting_invwarp.nii.gz
    # Strip point data that messes lmshoot up
    cat $LM_PROCRUSTES \
    | awk '$1 == "POINT_DATA" {stop=1} { if (stop!=1) print $0  }' \
    > $WITER/shooting_target_procrustes_nopd.vtk

    if [[ ! -f $SHOOTING_INVWARP ]]; then
    #Compute the invsere shooting
    time lmshoot -d 3 \
    -m $WITER/shooting_target_procrustes_nopd.vtk $LANDMARKS \
    -o $MOMENTA_INVERSE \
    -s 2.0 -l 5000 -n 40 -i 240 0 -f -O $WITER/movie/movie%04d.vtk

   # Convert the shooting result into a warp
   lmtowarp -d 3 -n 40 -r $REFSPACE \
   -m $MOMENTA_INVERSE -o $SHOOTING_INVWARP \
   -s 2.0
   fi
 fi

  # Warp the native space image into the template
  #reslice_subj $ID $REFSPACE $RESLICE_PTRN $SHOOTING_WARP $LM_PROCRUSTES_MAT,-1
  #reslice_subj_segsonly $ID $REFSPACE $RESLICE_PTRN $SHOOTING_WARP $LM_PROCRUSTES_MAT,-1

}


function fn_shooting_shape_avg()
{
  # The iteration of the loop
  iter=${1?}

  # The source landmarks
  SRC_LANDMARKS=${2?}

  # The mode - NCC or MST
  MODE=${3?}

  # The work dir
  if [[ $MODE == 'MST' ]]; then
    SHOOT_DIR=$WDIR/gshoot
  elif [[ $MODE == 'NCC' ]]; then
    SHOOT_DIR=$WDIR/gshoot_ncc
  fi

  WORK=$SHOOT_DIR/shape_avg/iter_${iter}

  # The reference space for shooting
  REFSPACE=$WDIR/template_mst/template_mst_csmed.nii.gz

  # The result landmarks - after shape averating
  SHAVG_LANDMARKS_NOPROC=$WORK/shavg_landmarks_noprocrustes.vtk
  SHAVG_LANDMARKS=$WORK/shavg_landmarks.vtk

  # Create the work dir
  mkdir -p $WORK

  # Average the momentum maps from the previous iteration
  /project/hippogang_2/pauly/bin/avgmesharr \
    $SHOOT_DIR/*/iter_${iter}/shooting_momenta.vtk \
    InitialMomentum $SRC_LANDMARKS $WORK/average_momenta.vtk

  # Perform the shooting and generate warp
  lmtowarp \
    -d 3 -n 40 -r $REFSPACE \
    -m $WORK/average_momenta.vtk \
    -o $WORK/average_momenta_warp.nii.gz -s 2.0

  # Apply the warp to the landmarks to bring them to shape-averaging position
  $GREEDY_HOME/greedy \
    -d 3 -rs $SRC_LANDMARKS $SHAVG_LANDMARKS_NOPROC \
    -rf $REFSPACE -r $WORK/average_momenta_warp.nii.gz

  # This transformation of the landmarks can cause shrinkage of the template. This
  # is not at all what we want in the template, we actually want the template to
  # keep its size during this iterative process. The way to correct this is to perform
  # procrustes between the source lanmarks and the new shape average
  vtkprocrustes $SRC_LANDMARKS $SHAVG_LANDMARKS_NOPROC $WORK/residual_procrustes.mat \
    | grep RMS_ | tee $WORK/procrustes_metric.txt
 # Applying the inverse of this procrustes to the SHAVG_LANDMARKS_NOPROC gives the new
  # template landmarks which are shape averaged but still the same size as the original
  # template.
  $GREEDY_HOME/greedy \
    -d 3 -rs $SRC_LANDMARKS $SHAVG_LANDMARKS \
    -rf $REFSPACE \
    -r $WORK/average_momenta_warp.nii.gz $WORK/residual_procrustes.mat,-1
}


function fn_paths_shooting_average()
{
  # Create the mean image and mean DB images
  mkdir -p $WDIR/template_gshoot
  for what in img csmed cslat phg srlm; do
    c3d $WDIR/gshoot/*/iter_$((GSHOOT_NITER-1))/reslice_*${what}.nii.gz \
      -mean -o $WDIR/template_gshoot/template_gshoot_${what}.nii.gz
  done
}

function fn_pair_gshoot_overlap()
{
  ID1=${1?}
  ID2=${2?}

  ITERDIR=iter_$((GSHOOT_NITER-1))

  WORK=$WDIR/template_gshoot/overlap/pairs
  mkdir -p $WORK

  # Compute Dice
  for what in phg srlm; do
    $C3D_HOME/c3d \
      $WDIR/gshoot/${ID1}/$ITERDIR/*_${what}.nii.gz \
      $WDIR/gshoot/${ID2}/$ITERDIR/*_${what}.nii.gz \
      -overlap 1 > $WORK/overlap_${what}_${ID1}_${ID2}.txt
  done

  # Compute NCC
  c3d \
    $WDIR/gshoot/${ID1}/$ITERDIR/*_img*.nii.gz \
    $WDIR/gshoot/${ID2}/$ITERDIR/*_img*.nii.gz \
    -ncc 4x4x4 -replace -NaN 0 NaN 0 -Inf 0 Inf 0 \
    $WDIR/template_gshoot/template_gshoot_phg.nii.gz -thresh 0.5 inf 1 0 \
    -lstat | awk 'NR==3 {print $2}' \
    > $WORK/ncc_${ID1}_${ID2}.txt
}

function fn_shooting_correction_iterative()
{
  MODE=${1?}

  START_ITER=0

  # The shooting dir
  if [[ $MODE == 'MST' ]]; then
    SHOOT_DIR=$WDIR/gshoot
    END_ITER=$GSHOOT_NITER
  elif [[ $MODE == 'NCC' ]]; then
    SHOOT_DIR=$WDIR/gshoot_ncc
    END_ITER=$GSHOOT_NITER_NCC
  fi

  mkdir -p $SHOOT_DIR/dump

  # Initial landmarks - sampled from template segmentation
  INIT_LANDMARKS=$SHOOT_DIR/template/root_landmarks.vtk

  # Loop $START_ITER
  for ((iter=$START_ITER;iter<$END_ITER;iter++)); do

    # Path to the landmarks at this iteration
    if [[ $iter -eq 0 ]]; then

      # Compute the initial landmarks
      $ROOT/scripts/job_launcher.sh -o $SHOOT_DIR/dump -N "shootinit" \
            "$0" fn_shooting_init_landmarks $MODE

      LANDMARKS=$INIT_LANDMARKS

    else

      LANDMARKS=$SHOOT_DIR/shape_avg/iter_$((iter-1))/shavg_landmarks.vtk

    fi

    # Launch all jobs
    N=$(cat $SUBJ_LIST | wc -l)
    for ((i=1;i<=$N;i++)); do
      $ROOT/scripts/job_launcher.sh -m 20G -o $SHOOT_DIR/dump -N "shoot_$i" \
	 "$0" fn_shooting_correction_subj_iter $i $iter $LANDMARKS $MODE
    done

    # Wait for completion
    $ROOT/scripts/job_launcher.sh -w "shoot_*" /bin/sleep 0

    # Launch the averaging job
    $ROOT/scripts/job_launcher.sh -m 30G -o $SHOOT_DIR/dump -N "avshoot_$i" \
          "$0" fn_shooting_shape_avg $iter $LANDMARKS $MODE

    # Wait for completion
    $ROOT/scripts/job_launcher.sh -w "avshoot_*" /bin/sleep 0


 done

}

function finalize_all()
{
  N=$(cat $SUBJ_LIST | wc -l)
  for ((i=1;i<=$N;i++)); do
     id=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
      $ROOT/scripts/job_launcher.sh -m 20G -o $WDIR/dump -N "final_${id}" "$0" finalize_subject $id
  done

  $ROOT/scripts/job_launcher.sh -w "final_*" /bin/sleep 0
}

# Compute forward and backward transforms to raw space and unwarp space for
# each subject. This is done relative to the last iteration of geodesic shooting
# i.e., the best template that we have been able to build.
function finalize_subject()
{
  ID=${1?}

  # Subject's preproc directory
  PDIR=$ROOT/preproc/${ID}

  # Subject's unwarped-space data
  read IMGUW CSMEDUW CSLATUW PHGUW SRLMUW <<<$(fn_input_unwarped_images $ID 1)

  # Subject's RAW MRI scan
  RAW_IMG=$PDIR/${ID}_raw_n4clip.nii.gz

  # Subject's source segmentation of DB and HF
  PHGSEG=$PDIR/${ID}_axisalign_phgsegshape_singlelabel.nii.gz
  SRLMSEG=$PDIR/${ID}_axisalign_srlmseg_sr.nii.gz

  # Subject's unwarp and dan/raw matrices
  UNWARP_MAT=$PDIR/unwarp/${ID}_unwarp.mat

<<"OLDCODE"
  # Extract format of axisalign to raw warpchain - saved in seg_to_raw.txt
  SEG_RAW_FORMAT=$(cat $ROOT/scripts/seg_to_raw.txt | awk -v id="$ID" '$1 == id {print $2}')

 if [[ $SEG_RAW_FORMAT -eq 1 ]]; then
        local SEG_TO_RAW_MAT="$PDIR/${ID}_warp.mat $PDIR/${ID}_inverse_axisalign.mat"
 fi

 if [[ $SEG_RAW_FORMAT -eq 2 ]]; then
        local SEG_TO_RAW_MAT="$PDIR/${ID}_inverse_axisalign.mat"
 fi
OLDCODE

  local SEG_TO_RAW_MAT="$PDIR/${ID}_transform_to_axisalign.mat,-1"

  # Subject's geodesic shooting transform from unwarp space to MST template space
  GSDIR=$WDIR/gshoot/${ID}/iter_$((GSHOOT_NITER-1))
  GSWARP=$GSDIR/shooting_warp.nii.gz
  GSWARP_INV=$GSDIR/shooting_invwarp.nii.gz
  GSPROC=$GSDIR/target_to_root_procrustes.mat

  # Subject's warp from shape-normalized space to NCC template - with the multilabel dataset, the NCC warp goes from template space to subject space
  # (Subject image was the fixed image)
  NCCDIR=$WDIR/template_ncc/iter$((NCC_NITER-1))
  NCCWARP_ROOT_INV=$NCCDIR/warp_root_template_to_${ID}.nii.gz
  NCCWARP_INV="$NCCWARP_ROOT_INV,64"

  ## Computing the inverse root was mesing up greedy during template to unwarp reslice (requested region error)
  #NCCWARP_INV="$NCCDIR/warp_template_to_${ID}.nii.gz"

  # Directory for the final output
  WORK=$WDIR/final/${ID}

  # Compute the warp using the root from subject to template space. Don't use root. Getting greedy error
  NCCWARP="$NCCWARP_ROOT_INV,-64"
  
<<"FIN"
  #This didn;t work well - inverse doesnt match up perfectly. Instead, computing lmshoot in the opposite direction
 GSWARP_INV=$WORK/${ID}_gshoot_invwarp_lmtowarp.nii.gz
  python VTKInverseMomenta.py -i $GSDIR/shooting_momenta.vtk -r $GSDIR/shooting_target_procrustes.vtk -o $WORK/${ID}_shooting_momenta_inv.vtk
  
  lmtowarp -d 3 -n 40 -r $WDIR/gshoot/${ID}/refspace.nii.gz \
  -m $WORK/${ID}_shooting_momenta_inv.vtk -o $GSWARP_INV \
  -s 2.0
FIN

# Warp chain from unwarped image space to template space
  CHAIN_UNWARPED_TO_TEMPLATE="$NCCWARP $GSWARP $GSPROC,-1"
  CHAIN_TEMPLATE_TO_UNWARPED="$GSPROC $GSWARP_INV $NCCWARP_INV" 

  # Warp chain from raw image space to template space
  CHAIN_RAW_TO_TEMPLATE="$CHAIN_UNWARPED_TO_TEMPLATE $UNWARP_MAT"

  # Warp chain from input image space (where segmentation drawn) to template
  CHAIN_INPUT_TO_TEMPLATE="$CHAIN_RAW_TO_TEMPLATE $SEG_TO_RAW_MAT"

  # Make the directory
  mkdir -p $WORK

  # Link the forward warps and other important data for completeness
  #ln -sf $NCCWARP_INV ${WORK}/${ID}_ncc_warp_inverse.nii.gz - the composite warp hasn't been generated as an image
  ln -sf $GSWARP ${WORK}/${ID}_gshoot_warp.nii.gz
  ln -sf $IMGUW ${WORK}/${ID}_unwarp_img.nii.gz
  ln -sf $NCCDIR/template_img.nii.gz $WORK/template_ncc.nii.gz

  # Write the chains to file
  echo $CHAIN_UNWARPED_TO_TEMPLATE > $WORK/chain_unwarped_to_template.txt
  echo $CHAIN_TEMPLATE_TO_UNWARPED > $WORK/chain_template_to_unwarped.txt
  echo $CHAIN_RAW_TO_TEMPLATE > $WORK/chain_raw_to_template.txt
  echo $CHAIN_INPUT_TO_TEMPLATE > $WORK/chain_input_to_template.txt

 #This doesn't work
  #$GREEDY_HOME/greedy -d 3 \
  #  -iw $GSWARP $GSWARP_INV -exp 4 -wp 0.0001
<<"DONE"

 # Create the best possible resampled intensity image (one interpolation less)
  $GREEDY_HOME/greedy -d 3 \
    -rf $WORK/template_ncc.nii.gz \
    -rm $RAW_IMG $WORK/${ID}_reslice_rawimg_to_template.nii.gz \
    -r $CHAIN_RAW_TO_TEMPLATE

 # Create the best possible HF and DB masks
  $GREEDY_HOME/greedy -d 3 \
    -ri LABEL 0.24vox \
    -rf $WORK/template_ncc.nii.gz \
    -rm $PHGSEG $WORK/${ID}_reslice_phgseg_to_template.nii.gz \
    -rm $SRLMSEG $WORK/${ID}_reslice_srlmseg_to_template.nii.gz \
    -r $CHAIN_INPUT_TO_TEMPLATE
DONE

}

function fn_template_ncc_loop
{
  # Which template are we building (ncc or aff_ncc)
  mode=${1?}

  # List of subjects
  N=$(cat $SUBJ_LIST | wc -l)

  # Iterations - number, starting, should be 0
  SITER=0

    # Perform the loop
    for ((iter=$SITER; iter<$NCC_NITER; iter++));do

    # Create a working directory for this iteration
    WORK=$WDIR/template_${mode}/iter${iter}
    mkdir -p $WORK/dump

    # Compute current average for each modality
    for mod in img phg csmed cslat srlm; do
		$ROOT/scripts/job_launcher.sh -o $WORK/dump -m 20G -N "average_${iter}_${mod}" \
        	"$0" fn_template_ncc_average $iter $mod $mode
    done

    # Wait for completion
    $ROOT/scripts/job_launcher.sh -w "average_${stage}*" /bin/sleep 0

    # register all images to the average and place in new folder for iter+1:
   for ((i=1;i<=$N;i++)); do
      	
	# submit ANTS job:
	$ROOT/scripts/job_launcher.sh -m 30G -o $WORK/dump -N "norm_${iter}_${i}" \
        "$0" fn_template_ncc_register ${i} ${iter} ${mode}
    	
	done

    # Wait for completion
        $ROOT/scripts/job_launcher.sh -w "norm_${stage}*" /bin/sleep 0
   done

}

# average images from iteration iter
function fn_template_ncc_average()
{

  # Iteration
  iter=${1?}

  # Modality
  mod=${2?}

  # Which template are we building (ncc or aff_ncc)
  mode=${3?}

  # Directory where to average
  WORK=$WDIR/template_${mode}/iter${iter}
  mkdir -p $WORK

  # List images to average
  if [[ $iter -eq 0 ]]; then

    if [[ $mode == "ncc" ]]; then
      # The inputs are in the paths directories
      INPUTS=$(ls $WDIR/gshoot/*/iter_$((GSHOOT_NITER-1))/reslice*${mod}*.nii.gz)
    elif [[ $mode == "aff_ncc" ]]; then
      INPUTS=$(ls $WDIR/template_${mode}/input/*/reslice*${mod}*.nii.gz)
    fi

  else

    # The inputs are outputs from last iteration
    INPUTS=$(ls $WDIR/template_${mode}/iter$((iter-1))/reslice*${mod}*.nii.gz)

  fi

  # Type of averaging to perform
  AVGMODE=0
  # This seems to make things worse
  if [[ $mod == "img" ]]; then AVGMODE=1; fi

  # Do the averaging
  AverageImages 3 $WORK/template_${mod}.nii.gz $AVGMODE ${INPUTS}

  # Scale template by 100 - why?
  if [[ $mod == "img" ]]; then
    c3d $WORK/template_${mod}.nii.gz -scale 100 -o $WORK/template_${mod}.nii.gz
  fi

  # Compute the mask
  if [[ $mod == "phg" ]]; then
    c3d $WORK/template_phg.nii.gz -thresh 0.5 inf 1 0 -dilate 1 10x10x10vox -o $WORK/template_mask.nii.gz
  fi
}

# This is the registration part of the groupwise intensity registration script.
# It can be run on full MST/GS output or on the affine part of it, the latter
# being a simulation of what conventional template building would produce on
# this dataset
function fn_template_ncc_register()
{

  # The index of the subject being registered using the MST
  idx=${1?}
  iter=${2?}

  # The mode for this function (ncc or aff_ncc)
  mode=${3?}

  # The work directory
  WORK=$WDIR/template_${mode}/iter${iter}

  # The ID of the subject
  ID=$(cat $SUBJ_LIST | head -n $idx | tail -n 1)

  # The native space images and registrations
  SEG_SHAPEONLY=0
  #read IMG_NAT CSMED_NAT CSLAT_NAT PHG_NAT SRLM_NAT <<<$(fn_input_unwarped_images $ID $SEG_SHAPEONLY)

  # The iteration path for gshoot
  GSDIR=$WDIR/gshoot/${ID}/iter_$((GSHOOT_NITER-1))

  # Build up the warp chain - all the 'shape' warps
  if [[ $mode == "ncc" ]]; then
    WARPCHAIN="$GSDIR/shooting_warp.nii.gz $GSDIR/target_to_root_procrustes.mat,-1"
    IMGMOV=$(ls $GSDIR/reslice_*img.nii.gz)
    PHGMOV=$(ls $GSDIR/reslice_*phg_int.nii.gz)
  elif [[ $mode == "aff_ncc" ]]; then
    WARPCHAIN="$WDIR/template_${mode}/input/$ID/proc_affine_${ID}.mat"
    IMGMOV=$(ls $WDIR/template_${mode}/input/$ID/reslice_*img.nii.gz)
    PHGMOV=$(ls $WDIR/template_${mode}/input/$ID/reslice_*phg_int.nii.gz)
  fi

  # Template
  TEMPLATE_IMG=$WORK/template_img.nii.gz

  # Pattern for reslicing
  RESLICE_PTRN=$WORK/reslice_to_template_${ID}_XXX.nii.gz

  # Create the working directory
  mkdir -p $WORK

 # Symlink the moving image in this directory
  ln -sf $IMGMOV $WORK/moving_${ID}_img.nii.gz

  # The output warp
  WARP=$WORK/warp_template_to_${ID}.nii.gz

  # The output warp - root
  WARP_ROOT=$WORK/warp_root_template_to_${ID}.nii.gz

  # Compute the mask over the moving image - only need to do this once
  PHGMOV_MASK=$WDIR/template_${mode}/moving_${ID}_mask.nii.gz
 if [[ $iter -eq 0 ]]; then
  c3d $PHGMOV -thresh 1 1 1 0 -dilate 1 10x10x10vox -as VISIBLE $PHGMOV -thresh 2 2 1 0 -dilate 1 1x1x1vox -insert VISIBLE  2 -scale -1 -add -o $PHGMOV_MASK
fi

  # Perform registration with NCC -
  # Greedy only supports a mask over the fixed image. Reverse the order of the registration
  # so that the template is warped to the subject space. The subject can then be warped to the template
  # space using the inverse warp ( root,-64)

  #Should be NCC 2x2x2 - why is it crashing??
 $GREEDY_HOME/greedy -d 3 \
    -w 1 -i $IMGMOV $TEMPLATE_IMG \
    -gm $PHGMOV_MASK \
    -m NCC 2x2x2 \
    -n 100x100x50 -s 0.6mm 0.2mm -e 0.5 \
    -threads $NSLOTS -oroot $WARP_ROOT \
    -sv -wp 0

# to save storage, don't need to save full warp here. Only the root.

<< "TEMPLATE_FIXED"
 $GREEDY_HOME/greedy -d 3 \
    -w 1 -m NCC 2x2x2 -i $TEMPLATE_IMG $IMGMOV \
    -gm $WORK/template_mask.nii.gz \
    -n 100x100x50 -s 0.6mm 0.2mm -e 0.5 \
    -threads $NSLOTS \
    -sv -wp 0 -oroot $WARP_ROOT \
    -o $WARP
TEMPLATE_FIXED

  # Apply registration at this iteration
  reslice_subj $ID $TEMPLATE_IMG $RESLICE_PTRN "$WARP_ROOT,-64" $WARPCHAIN
}

function fn_template_ncc_all_pairs_overlap
{

  iter=${1?}
  mode=${2?}

  N=$(cat $SUBJ_LIST | wc -l)
  for ((i=1;i<=$N;i++)); do
    IDi=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
    for ((j=$((i+1));j<=$N;j++)); do
      if [[ $i -ne $j ]]; then
        IDj=$(cat $SUBJ_LIST | head -n $j | tail -n 1)
	$ROOT/scripts/job_launcher.sh -m 20G -o $WDIR/dump -N "ovl_${IDi}_${IDj}" \
          "$0" fn_template_ncc_pair_overlap $IDi $IDj $iter $mode
      fi
    done
  done

  # Wait for completion
  $ROOT/scripts/job_launcher.sh -w "ovl_*" /bin/sleep 0

  $ROOT/scripts/job_launcher.sh -m 40G -o $WDIR/dump -N "avg_ncc" \
          "$0" fn_template_ncc_average_nccmap

}

function fn_template_ncc_pair_overlap()
{

  ID1=${1?}
  ID2=${2?}
  iter=${3?}
  mode=${4?}

  WORK=$WDIR/template_${mode}/iter${iter}
  TEMPLATE_MASK=$WDIR/template_${mode}/iter${iter}/template_mask.nii.gz

  mkdir -p $WORK/overlap/pairs

  # Compute Dice
  for what in phg srlm; do
    $C3D_HOME/c3d \
      $WORK/reslice_to_template_${ID1}_${what}.nii.gz \
      $WORK/reslice_to_template_${ID2}_${what}.nii.gz \
      -overlap 1 > $WORK/overlap/pairs/overlap_${what}_${ID1}_${ID2}.txt
  done

  # Compute NCC
  c3d \
    $WORK/reslice_to_template_${ID1}_img.nii.gz \
    $WORK/reslice_to_template_${ID2}_img.nii.gz \
    -ncc 4x4x4 -replace -NaN 0 NaN 0 -Inf 0 Inf 0 \
    $WORK/template_phg.nii.gz -thresh 0.5 inf 1 0 \
    -lstat | awk 'NR==3 {print $2}' \
    > $WORK/overlap/pairs/ncc_${ID1}_${ID2}.txt

    $GREEDY_HOME/greedy -d 3 \
    -i $WORK/reslice_to_template_${ID1}_img.nii.gz $WORK/reslice_to_template_${ID2}_img.nii.gz \
    -gm $TEMPLATE_MASK -metric -m NCC 4x4x4 -o \
   $WORK/overlap/pairs/ncc_map_${ID1}_${ID2}.nii.gz
   
}


fn_template_ncc_average_nccmap()
{

   WORK=$WDIR/template_ncc/iter5
   $GREEDY_HOME/greedy_template_average -d 3 \
   -i $WORK/overlap/pairs/ncc_map_*.nii.gz $WORK/overlap/ncc_map.nii.gz

}

function misc_final_stats()
{

<<"DONE"
  local TEMPLATE_PATHS="
    /data/jux/sravikumar/atlasPHG2019/mst_multilabel_v2/template_mst/overlap/pairs \
    /data/jux/sravikumar/atlasPHG2019/mst_multilabel_v2/template_gshoot/overlap/pairs \
    /data/jux/sravikumar/atlasPHG2019/mst_multilabel_v2/template_ncc/iter5/overlap/pairs \
   /data/jux/sravikumar/atlasPHG2019/mst_multilabel_v2/template_aff_ncc/iter5/overlap/pairs"
for P in $TEMPLATE_PATHS; do
    summarize_template_stats $P
  done
DONE

# Launch all jobs
  for mode in mst ncc gshoot aff_ncc; do
    qsubp2 -cwd -o $WDIR/dump -j y -N "slice_${mode}" \
          $0 template_screenshots $mode
  done
  # Wait for completion
  qsub -cwd -o $WDIR/dump -j y -hold_jid "slice_*" -sync y -b y sleep 1
}


function summarize_template_stats()
{
  WORK=${1?}

  cat $WORK/overlap_phg_* | awk -F ',' '{print $5}' > $TMPDIR/phg.txt
  cat $WORK/overlap_srlm_* | awk -F ',' '{print $5}' > $TMPDIR/srlm.txt
  cat $WORK/ncc_*.txt > $TMPDIR/ncc.txt

  echo "OVL_PHG OVL_SRLM NCC" > $WORK/pairwise_stats.txt
  paste $TMPDIR/phg.txt  $TMPDIR/srlm.txt $TMPDIR/ncc.txt >> $WORK/pairwise_stats.txt
}

function template_screenshots()
{

  local mode=${1?}
  local WORK=$WDIR/final/templates/template_${mode}
  local CROPSPACE=$ROOT/scripts/manual/cropspace.nii.gz

  mkdir -p $WORK/reslice $WORK/slices/template/

  # Build up the list of images to reslice - use links
  for id in $(cat $SUBJ_LIST); do

    local RESLICE
    local RESLICE_PHG
    case "$mode" in
      ncc)
        RESLICE=$WDIR/template_ncc/iter$((NCC_NITER-1))/reslice_to_template_${id}_img.nii.gz
        ;;
      aff_ncc)
        RESLICE=$WDIR/template_aff_ncc/iter$((NCC_NITER-1))/reslice_to_template_${id}_img.nii.gz
        ;;
      aff_only)
        RESLICE=$WDIR/template_aff_ncc/iter0/moving_${id}_img.nii.gz
        ;;
      mst)
        RESLICE=$(ls $WDIR/paths/${id}/final/reslice_${id}_*img.nii.gz)
	RESLICE_PHG=$(ls $WDIR/paths/${id}/final/reslice_${id}_*phg.nii.gz)
        ;;
gshoot)
        RESLICE=$(ls $WDIR/gshoot/${id}/iter_$((GSHOOT_NITER-1))/reslice*img.nii.gz)
        ;;
    esac

    ln -sf $RESLICE $WORK/reslice/reslice_${id}.nii.gz
    ln -sf $RESLICE_PHG $WORK/reslice/reslice_phg_${id}.nii.gz
  done

  if [[ ! -f $WORK/template_${mode}.nii.gz ]]; then
  # Create templates from the images by averaging
  AverageImages 3 $WORK/template_${mode}.nii.gz 0 $WORK/reslice/*.nii.gz
  fi

  c3d \
    $CROPSPACE $WORK/template_${mode}.nii.gz \
    -reslice-identity -o $WORK/crop_template_${mode}.nii.gz

  # Generate slices for the template
  c3d \
    $WORK/template_${mode}.nii.gz \
    -stretch 0 850 0 255 -clip 0 255 -popas X \
    -type uchar \
    -push X -slice x 376 -flip x -flip y -o $WORK/slices/coronal_ant_${mode}_template.png \
    -push X -slice y 126  -flip x -flip y -o $WORK/slices/sagittal_${mode}_template.png \
    -push X -slice z 108  -flip x -flip y -o $WORK/slices/axial_${mode}_template.png
  
   c3d \
    $WORK/template_${mode}.nii.gz \
    -stretch 0 850 0 255 -clip 0 255 -popas X \
    -type uchar \
    -push X -slice x 387 -flip x -flip y -o $WORK/slices/coronal_1_${mode}_template.png \
    -push X -slice x 327 -flip x -flip y -o $WORK/slices/coronal_2_${mode}_template.png \
    -push X -slice x 321 -flip x -flip y -o $WORK/slices/coronal_3_${mode}_template.png \
    -push X -slice x 292 -flip x -flip y -o $WORK/slices/coronal_4_${mode}_template.png \
    -push X -slice x 260 -flip x -flip y -o $WORK/slices/coronal_5_${mode}_template.png \
    -push X -slice x 213 -flip x -flip y -o $WORK/slices/coronal_6_${mode}_template.png

  # Generate slices for individual images y/z =  #138/128 for ncc gshoot
  for id in $(cat $SUBJ_LIST); do

 #coronal was 314
    c3d \
      $WORK/reslice/reslice_${id}.nii.gz \
      -stretch 0 850 0 255 -clip 0 255 -popas X \
      -type uchar \
      -push X -slice x 376 -flip x -flip y -o $WORK/slices/coronal_ant_${mode}_${id}.png \
      -push X -slice y 126  -flip x -flip y -o $WORK/slices/sagittal_${mode}_${id}.png \
      -push X -slice z 108  -flip x -flip y -o $WORK/slices/axial_${mode}_${id}.png

<<"NONEED"
    c3d \
       $WORK/reslice/reslice_${id}.nii.gz \
      -stretch 0 850 0 255 -clip 0 255 -popas X \
      -type uchar \
      -push X -slice x 387 -flip x -flip y -o $WORK/slices/coronal_1_${mode}_${id}.png \
      -push X -slice x 327 -flip x -flip y -o $WORK/slices/coronal_2_${mode}_${id}.png \
      -push X -slice x 321 -flip x -flip y -o $WORK/slices/coronal_3_${mode}_${id}.png \
      -push X -slice x 292 -flip x -flip y -o $WORK/slices/coronal_4_${mode}_${id}.png \
      -push X -slice x 260 -flip x -flip y -o $WORK/slices/coronal_5_${mode}_${id}.png \
      -push X -slice x 213 -flip x -flip y -o $WORK/slices/coronal_6_${mode}_${id}.png 

     c3d \
       $WORK/reslice/reslice_phg_${id}.nii.gz \
      -stretch 0 1 0 255 -clip 0 255 -popas X \
      -type uchar \
      -push X -slice x 387 -flip x -flip y -o $WORK/slices/coronal_1_${mode}_${id}_phg.png \
      -push X -slice x 327 -flip x -flip y -o $WORK/slices/coronal_2_${mode}_${id}_phg.png \
      -push X -slice x 327 -flip x -flip y -o $WORK/slices/coronal_2_${mode}_${id}_phg.png \
      -push X -slice x 321 -flip x -flip y -o $WORK/slices/coronal_3_${mode}_${id}_phg.png \
      -push X -slice x 292 -flip x -flip y -o $WORK/slices/coronal_4_${mode}_${id}_phg.png \
      -push X -slice x 260 -flip x -flip y -o $WORK/slices/coronal_5_${mode}_${id}_phg.png \
      -push X -slice x 213 -flip x -flip y -o $WORK/slices/coronal_6_${mode}_${id}_phg.png
NONEED

done


}


function fn_prepare_affine_template_ncc_all()
{
  # Launch all jobs
  N=$(cat $SUBJ_LIST | wc -l)
  for ((i=1;i<=$N;i++)); do
    qsubp4 -cwd -o $WDIR/dump -j y -N "affprep_$i" \
          $0 fn_prepare_affine_template_ncc $i
  done

  # Wait for completion
  qsub -cwd -o $WDIR/dump -j y -hold_jid "affprep_*" -sync y -b y sleep 1
}


# Determine if the subject with given index needs flipping when doing procrustes
# alignment with the MST root
function fn_flip_to_root()
{
  local idx=${1?}
  local ID=$(cat $ROOT/scripts/subj24_08072019.txt | head -n $idx | tail -n 1)
  local SIDE=$(echo $ID | awk -F '-' '{print $NF}')
  local flip=""
  OLD_ATLAS_SIDE="L"

  if [[ $SIDE != $OLD_ATLAS_SIDE ]]; then flip="-f"; fi

  echo $flip
}


# Prepare for affine-then-intensity template population build
function fn_prepare_affine_template_ncc()
{

 export LD_LIBRARY_PATH=/data/picsl/pauly/lib:/data/jux/sravikumar/packages/vtk-release-6.3/bin/lib:/data/jux/sravikumar/packages/glibc-2.14/lib:/data/jux/sravikumar/atlasPHG2019/pkgs/libc:$LD_LIBRARY_PATH

  idx=${1?}

  # The ID of the subject
  ID=$(cat $SUBJ_LIST | head -n $idx | tail -n 1)

  SEG_SHAPEONLY=1
  # The native space images and registrations
  read IMGNAT CSMEDNAT CSLATNAT PHGNAT SRLMNAT <<<$(fn_input_unwarped_images $ID $SEG_SHAPEONLY)

  # The iteration path for gshoot
  GSDIR=$WDIR/gshoot/${ID}/iter_$((GSHOOT_NITER-1))

  # The fixed and moving landmarks for gshoot
  GSL_FIXED=$GSDIR/landmarks.vtk
  GSL_MOVING=$GSDIR/../shooting_target_native.vtk

  # Do we need a flip?
  flip=$(fn_flip_to_root $idx)

  # Directory for the work
  WORK=$WDIR/template_aff_ncc/input/${ID}

  # The affine matrix
  AFFPROC=$WORK/proc_affine_${ID}.mat

  # The reference image where to do registrations
  REFSPACE=$WDIR/template_gshoot/template_gshoot_phg.nii.gz

  mkdir -p $WORK

  # Extract the affine component of the geodesic shooting
  vtkprocrustes $flip $GSL_FIXED $GSL_MOVING $AFFPROC

  # Apply transformation to the images
  reslice_subj $ID $REFSPACE $WORK/reslice_to_template_${ID}_XXX.nii.gz $AFFPROC
  reslice_subj_segsonly $ID $REFSPACE $WORK/reslice_to_template_${ID}_XXX.nii.gz $AFFPROC

}


function main()
{
  
  mkdir -p $WDIR/dump

  # Compute the MST
  fn_all_pairwise

  fn_print_adjacency

  # MST root is already manually axisaligned.
  #Instead pad root image and warp to same space as axisaligned root image
  $ROOT/scripts/job_launcher.sh -m 15G -o $WDIR/dump -N "root_to_oldatlas" "$0" fn_register_root_to_oldatlas

  # Perform registration along MST paths and compute template and template-space overlaps
  fn_all_paths fn_register_mst
  $ROOT/scripts/job_launcher.sh -m 40G -o $WDIR/dump -N "avg_path" "$0" fn_paths_average
  
  #Generate screenshots here instead of in the final stats step to check the atlas quality
  template_screenshots mst
  fn_all_pair_paths_overlap fn_pair_paths_overlap

  # Apply the geodesic shooting correction and reevaluate the overlaps and such
  fn_shooting_correction_iterative MST
  
  $ROOT/scripts/job_launcher.sh -m 40G -o $WDIR/dump -N  "avg_path" " $0" fn_paths_shooting_average
  fn_all_pair_paths_overlap fn_pair_gshoot_overlap

  # Check to see if all the subjects look ok
  template_screenshots gshoot

 # Build a template in intensity space - consecutive averaging and registration
  fn_template_ncc_loop ncc

 fn_template_ncc_all_pairs_overlap $((NCC_NITER-1)) ncc

 fn_shooting_correction_iterative NCC

 finalize_all

# Build a comparison template using intensity groupwise registration after only
# affine alignment - no MST/GS
fn_prepare_affine_template_ncc_all
fn_template_ncc_loop aff_ncc
fn_template_ncc_all_pairs_overlap $((NCC_NITER-1)) aff_ncc

# Compute all the summary stats
qsubp4 -cwd -o $WDIR/dump -j y -N "misc_stats" -sync y $0 misc_final_stats

}


if [[ $1 ]]; then

  command=$1
  shift
        $command $@

else

  main

fi
