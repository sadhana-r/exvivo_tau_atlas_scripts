#!/bin/bash
set -x -e

ROOT=/project/hippogang_3/sravikumar/atlasPHG2019
source $ROOT/scripts/common.sh


# Relies on list of subjects manually entered in manifest
function copy_hires_mri()
{
  while IFS=$'\t' read -r ID_SR ID FWPATH args; do

    echo $ID_SR
    # Create the input directory
    IDIR=$ROOT/pmatlas_raw/$ID/hires_mri
    mkdir -p $IDIR

    # Go there
    pushd $IDIR > /dev/null

    # Base filename
    local FN=$(basename "$FWPATH")

    # Check for existing file
    if [[ -f $FN || -f ${FN}.gz ]]; then
      echo Skipping $ID:$FN
      continue
    fi

    # Copy needed files
    fw download -o $FN "$FWPATH"

    # Compress if needed
    if [[ $FN =~ nii$ ]]; then
      gzip $FN
      ln -sf ${FN}.gz ${ID}_mri_hires.nii.gz
    else
      ln -sf $FN ${ID}_mri_hires.nii.gz
    fi

    popd > /dev/null

  done < $ROOT/scripts/manual/hiresmri_manifest.txt
}

function map_to_axisalign()
{

 SEQ=($(seq -f "%04g" 0 50))

 #i=15
 while IFS=$'\t' read -r ID_SR ID FWPATH args; do
 #if grep -Fxq "$ID_SR" $ROOT/scripts/new_cases_july2022.txt
 #then
   # Create the input directory
   IDIR=$ROOT/pmatlas_raw/${ID}/hires_n4

   #Create the output directory
   ODIR=$ROOT/pmatlas_raw/$ID/hires_n4_axisalign
   mkdir -p $ODIR

   HIRES_MRI_RAW=${IDIR}/${ID_SR}_mri_raw_n4_clip.nii.gz
   TRANSFORM_TO_AXISALIGN=$ROOT/scripts/manual/raw_to_axisalign_transforms/${ID_SR}_raw_to_axisalign.mat
   HIRES_AXISALIGN=$ODIR/${ID_SR}_n4clip_axisalign_img.nii.gz
	
   REF_SPACE=$ROOT/atlas2016/template_img_pad.nii.gz

   if [[ ! -f $TRANSFORM_TO_AXISALIGN ]]; then

	  continue
   else

   	if [[ ! -f $HIRES_AXISALIGN  ]]; then   	
		greedy -d 3 \
   		-rf $REF_SPACE \
   		-rm $HIRES_MRI_RAW $HIRES_AXISALIGN \
   		-ri LINEAR \
   		-r $TRANSFORM_TO_AXISALIGN
   	fi
	
   fi
       #cp $HIRES_AXISALIGN $ODIR/${ID_SR}_${SEQ[$i]}_0000.nii.gz
       #c3d $HIRES_AXISALIGN -thresh 0 inf 0 0 -o $ODIR/${ID_SR}_${SEQ[$i]}_0001.nii.gz
       i=$((i+1))
   #fi
done < $ROOT/scripts/manual/hiresmri_manifest.txt

}

function n4correct()
{


  #SEQ=($(seq -f "%04g" 0 50))
  #i=0

  while IFS=$'\t' read -r ID_SR ID FWPATH args; do
  
  #if [[ "$ID_SR" == "HNL29_18-L" ]];
  #if grep -Fxq "$ID_SR" $ROOT/scripts/subj_aaic2022.txt
  #then

  IDIR=$ROOT/pmatlas_raw/$ID/hires_mri
 
  #Create the output directory
  #NNUNET_DIR=$ROOT/hires_fornnunet
  #mkdir -p $NNUNET_DIR
  ODIR=$ROOT/pmatlas_raw/$ID/hires_n4
  mkdir -p $ODIR

  HIRES_MRI_RAW=${IDIR}/${ID}_mri_hires.nii.gz
  PRECLIP=/tmp/preclip.nii.gz
  HIRES_N4=$ODIR/${ID_SR}_mri_raw_n4.nii.gz
  HIRES_N4_CLIP=$ODIR/${ID_SR}_mri_raw_n4_clip.nii.gz  
  
  if [[ ! -f $HIRES_N4_CLIP ]]; then
#	
#	echo "Done"
#	cp $HIRES_N4_CLIP $NNUNET_DIR/${ID_SR}_${SEQ[$i]}_0000.nii.gz
#        i=$((i+1))
#	continue
#  fi

  # Perform N4
  $C3D_HOME/c3d $HIRES_MRI_RAW -stretch 0.1% 99.9% 0 1000 -clip 0 1000 -o $PRECLIP
  
  N4BiasFieldCorrection -d 3 -i $PRECLIP -o $HIRES_N4

  # Rescale the raw image to 0 - 1000 range (as input images already have)
  $C3D_HOME/c3d $HIRES_N4 -stretch 0.1% 99.9% 0 1000 -clip 0 1000 -o $HIRES_N4_CLIP

  else
	echo "N4 corrected image exists"

  fi
  #fi

  done < $ROOT/scripts/manual/hiresmri_manifest.txt
}


#copy_hires_mri
#n4correct

#Before running this next step, have to manually generate the *_raw_to_axisalign.mat transform. This is the ssaved to
# ../manual/raw_to_axisalign_transforms. This transform is the rigid trasnformation between the subject's n4_clip scan and
# template_img_pad.nii.gz which is also saved in ../manual
map_to_axisalign
