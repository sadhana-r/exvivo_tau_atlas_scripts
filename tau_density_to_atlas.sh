#!/bin/bash
set -x -e

ATLAS_DIR=/project/hippogang_3/sravikumar/atlasPHG2019

source $ATLAS_DIR/scripts/common.sh

TAU_DIR=/project/hippogang_2/pauly/tau_atlas

WORK_DIR=$ATLAS_DIR/tau_atlas
mkdir -p $WORK_DIR


function vis_to_raw_space()
{

   ID=${1?}
   ID_PY=$(cat $ATLAS_DIR/scripts/subj_tau_atlas_2022.txt | awk -v id=$ID '{if ($1 == id) print $2}')

   WDIR=$WORK_DIR/$ID

   #Input files and warps from Paul's folder
   TAU_WDIR=$TAU_DIR/work/$ID_PY
   #RAW_MRI=$(find $ATLAS_DIR/pmatlas_raw/ -name "*116748*_wrongheader.nii.gz")
   RAW_MRI=$ATLAS_DIR/preproc/$ID/${ID}_raw.nii.gz
   TAU_MAP=$TAU_WDIR/recon_native/raw_hires/${ID_PY}_density_times_softmask_Tau_tangles.nii.gz
   TAU_MAP_MASK=$TAU_WDIR/historeg/whole/${ID_PY}_mask_Tau_tangles.nii.gz
   VIS_MRI=$TAU_WDIR/mri/${ID_PY}_mri_hires_vis.nii.gz


   #Registration from raw to axisalign space
   RAW_TO_AXISALIGN_AFFINE=$ATLAS_DIR/preproc/${ID}/${ID}_transform_to_axisalign.mat
   AXISALIGN_IMG=$ATLAS_DIR//inputs/${ID}/${ID}_axisalign_img.nii.gz

   if [[ -f $TAU_MAP ]]; then

   mkdir -p $WDIR
   ln -sf $TAU_MAP $WDIR/${ID}_density_tau_tangles_vis.nii.gz
   ln -sf $VIS_MRI $WDIR/${ID}_mri_vis.nii.gz
   ln -sf $TAU_MAP_MASK $WDIR/${ID}_density_mask_tau_tangles_vis.nii.gz

   # Output files
   TAU_MAP_RAW=$WDIR/${ID}_density_tau_tangles_raw.nii.gz
   VIS_TO_RAW=$WDIR/${ID}_mri_vis_to_raw.nii.gz
   RAW_TO_VIS=$WDIR/${ID}_mri_raw_to_vis_test.nii.gz
   TAU_MAP_AXISALIGN=$WDIR/${ID}_density_tau_tangles_axisalign.nii.gz
   TAU_MAP_MASK_AXISALIGN=$WDIR/${ID}_density_mask_tau_tangles_axisalign.nii.gz
   TAU_MAP_SMOOTH_AXISALIGN=$WDIR/${ID}_density_tau_tangles_smooth_axisalign.nii.gz

    # This is a special case. the raw scan used in the histology recon has the wrong header info.
   # The corrected raw scan was used during atlas construction. Need to inlclude an additional regitration step
   # between these two raw scans. A smilar situation is likely for HNL11-15 when we come to it.
  if [[ "${ID_PY}" == "INDD116748" || "${ID_PY}" == "HNL-11-15" ]];then

   id=$( echo $ID_PY | cut -c 5-)
   HISTO_RAW=$(find $ATLAS_DIR/pmatlas_raw/ -name "*${id}*_wrongheader*.nii.gz")
   MAT=$WDIR/${ID}_mri_atlas_to_histo_raw.mat

    if [[ ! -f $MAT ]]; then
        greedy -d 3 -a \
        -i $HISTO_RAW $RAW_MRI \
        -dof 6 -n 100x50x0 -m NCC 4x4x4 \
        -o $MAT
    fi

   WARP_CHAIN_VIS_TO_RAW="$MAT,-1 $MOLD_TO_HIRES_WARP $HIRES_TO_MOLD_AFFINE,-1 $MOLD_REORIENT_VIS,-1"

 else

   WARP_CHAIN_VIS_TO_RAW="$MOLD_TO_HIRES_WARP $HIRES_TO_MOLD_AFFINE,-1 $MOLD_REORIENT_VIS,-1"

 fi

 echo $WARP_CHAIN_VIS_TO_RAW > $WDIR/chain_vis_to_raw.txt

 WARP_CHAIN_VIS_TO_AA="$RAW_TO_AXISALIGN_AFFINE $WARP_CHAIN_VIS_TO_RAW"
 echo $WARP_CHAIN_VIS_TO_AA > $WDIR/chain_vis_to_axisalign.txt

<<"SKIP"
  # Test run the forward warp
   $GREEDY_HOME/greedy -d 3 \
   -rf $VIS_MRI \
   -rm $RAW_MRI  $RAW_TO_VIS \
   -r $MOLD_TO_HIRES_INVWARP $HIRES_TO_MOLD_AFFINE $MOLD_REORIENT_VIS
SKIP

   $GREEDY_HOME/greedy -d 3 \
   -rf $RAW_MRI \
   -rm $TAU_MAP $TAU_MAP_RAW \
   -r $WARP_CHAIN_VIS_TO_RAW

   $GREEDY_HOME/greedy -d 3 \
   -rf $RAW_MRI \
   -rm $VIS_MRI $VIS_TO_RAW \
   -r $WARP_CHAIN_VIS_TO_RAW

<<"FIN"
    $GREEDY_HOME/greedy -d 3 \
   -rf $AXISALIGN_IMG \
   -rm $WDIR/${ID}_density_tau_tangles_vis.nii.gz $TAU_MAP_AXISALIGN \
   -rm $WDIR/${ID}_density_mask_tau_tangles_vis.nii.gz $TAU_MAP_MASK_AXISALIGN \
   -r $WARP_CHAIN_VIS_TO_AA

  c3d $TAU_MAP_MASK_AXISALIGN -clip 0 inf -shift 0.001 \
  -as M -smooth-fast 0.5mm -reciprocal -as RSM \
  $TAU_MAP_AXISALIGN -push M -times -smooth-fast 0.5mm \
  -times -clip 0 5 -info -o $TAU_MAP_SMOOTH_AXISALIGN
FIN
fi

}

warp_tau_to_template()
{

  export LD_LIBRARY_PATH=/data/picsl/pauly/lib:/data/jux/sravikumar/packages/vtk-release-6.3/bin/lib:/data/jux/sravikumar/packages/glibc-2.14/lib:/data/jux/sravikumar/atlasPHG2019/pkgs/libc:$LD_LIBRARY_PATH

  ID=${1?}
  SDIR=$ATLAS_DIR/mst_multilabel_v2/final/$ID

  ID_PY=$(cat $ATLAS_DIR/scripts/tau_atlas_subj.txt | awk -v id=$ID '{if ($1 == id) print $2}')
  TAU_WDIR=$TAU_DIR/work/$ID_PY
  TAU_MAP=$WORK_DIR/${ID}/${ID}_density_tau_tangles_vis.nii.gz
  TAU_MAP_MASK=$WORK_DIR/${ID}/${ID}_density_mask_tau_tangles_vis.nii.gz

  #Template location
  TEMPLATE_MRI=$ATLAS_DIR/mst_multilabel_v2/final/templates/template_ncc/template_ncc.nii.gz
  CHAIN_VIS_TO_RAW=$( cat $WORK_DIR/$ID/chain_vis_to_raw.txt)
  CHAIN_RAW_TO_TEMPLATE=$( cat $SDIR/chain_raw_to_template.txt)

  TAU_TEMPLATE=$WORK_DIR/${ID}/${ID}_template_density_tau_tangles.nii.gz
  TAU_TEMPLATE_MASK=$WORK_DIR/${ID}/${ID}_template_density_mask_tau_tangles.nii.gz
  TAU_TEMPLATE_SMOOTH=$WORK_DIR/${ID}/${ID}_template_density_smooth_tau_tangles.nii.gz

  #Map tau to template
  $GREEDY_HOME/greedy -d 3 \
  -rf $TEMPLATE_MRI \
  -rm $TAU_MAP $TAU_TEMPLATE \
  -rm $TAU_MAP_MASK $TAU_TEMPLATE_MASK \
  -r $CHAIN_RAW_TO_TEMPLATE $CHAIN_VIS_TO_RAW

  c3d $TAU_TEMPLATE_MASK -clip 0 inf -shift 0.001 \
  -as M -smooth-fast 0.5mm -reciprocal -as RSM \
  $TAU_TEMPLATE -push M -times -smooth-fast 0.5mm \
  -times -clip 0 5 -info -o $TAU_TEMPLATE_SMOOTH
}


function tau_map_average()
{

 export LD_LIBRARY_PATH=/data/picsl/pauly/lib:/data/jux/sravikumar/packages/vtk-release-6.3/bin/lib:/data/jux/sravikumar/packages/glibc-2.14/lib:/data/jux/sravikumar/atlasPHG2019/pkgs/libc:$LD_LIBRARY_PATH

 N=$(cat $ATLAS_DIR/scripts/tau_atlas_subj.txt | wc -l)

 MERGE_LIST=""

 for((i=10;i<=$N;i++)); do

   ID=$(cat $ATLAS_DIR/scripts/tau_atlas_subj.txt | head -n $i | tail -n 1 | awk '{print $1}')
   include=$(cat $ATLAS_DIR/scripts/tau_atlas_subj.txt | awk -v id=$ID '{if ($1 == id) print $3}')
   if [[ $include == "Y" ]]; then
      MERGE_ENTRY=$WORK_DIR/${ID}/${ID}_template_density_smooth_tau_tangles.nii.gz

      MERGE_LIST="$MERGE_LIST $MERGE_ENTRY"
   fi
 done

 ${ANTSDIR}/AverageImages 3 $WORK_DIR/template_average_tau_nofourR_smoothed.nii.gz 0 $MERGE_LIST

}

function analyze_density()
{
  TSEG=${1?}
  TLAB=${2?}
  OUTCSV=${3?}
  PERCENTILE=${4?}

  echo "Specimen,Label,LVox,LMaskVox,Tau" > $OUTCSV

  N=$(cat $ATLAS_DIR/scripts/subj_tau_atlas.txt | wc -l)

  for((i=1;i<=$N;i++)); do

    id=$(cat $ATLAS_DIR/scripts/subj_tau_atlas.txt | head -n $i | tail -n 1 | awk '{print $1}')

    DENSITY=$WORK_DIR/${id}/${id}_template_density_tau_tangles.nii.gz
    MASK=$WORK_DIR/${id}/${id}_template_density_mask_tau_tangles.nii.gz

    # Mask-normalized density map
    N_DENSITY=$WORK_DIR/${id}/${id}_template_density_tau_tangles_norm.nii.gz

    # Binary mask of where normalization can be trusted
    N_MASK=$WORK_DIR/${id}/${id}_template_density_mask_tau_tangles_bin_mask.nii.gz

<<"DONE"
    # Process to generate a normalized density mask
    $C3D_HOME/c3d $DENSITY $MASK \
      -foreach -smooth-fast 0.5mm -endfor -as M \
      -reciprocal -times \
      -push M -thresh 0.25 inf 1 0 -o $N_MASK \
      -times -replace NaN 0 Inf 0 -Inf 0 -o $N_DENSITY
DONE

  # Get a list of all qualifying labels - ones that contain at least 25% mask
    ALL_LABELS=$($C3D_HOME/c3d $N_MASK $TSEG -lstat | awk 'NR>1 && $1 > 0 && $2 > 0.25 {print $1}')

    # Process to generate 90th percentiles over anatomical labels
    for lab in $ALL_LABELS; do

      LABSTR=$(cat $TLAB | awk -v l=$lab '$1==l {print $8}' | awk -F'"' '{print $2}')

<<"SKIP"
      #Used unnormalized for unfolding work to be consistent
      $C3D_HOME/c3d $DENSITY -shift 1 \
        $N_MASK $TSEG \
        -thresh $lab $lab 1 0 -voxel-sum -times -voxel-sum -as M -times \
        -pim fq -clip 0 90% -push M -lstat > /tmp/out.txt
SKIP

       #Dilate hotspot seg a little
      $C3D_HOME/c3d $N_DENSITY -shift 1 \
      $N_MASK $TSEG \
      -thresh $lab $lab 1 0 -dilate 1 3x3x1 -voxel-sum -times -voxel-sum -as M -times \
      -pim fq -clip 0 90% -push M -lstat > /tmp/out.txt

      #Dilate hotspot seg a little. Mean across all density values (without 90% clipping)
      $C3D_HOME/c3d $N_DENSITY -shift 1 \
      $N_MASK $TSEG \
      -thresh $lab $lab 1 0 -dilate 1 3x3x1 -voxel-sum -times -voxel-sum -as M -times \
      -push M -lstat > /tmp/out_mean.txt

      NLAB=$(cat /tmp/out.txt | awk 'NR==1 {print $3}')
      NMSK=$(cat /tmp/out.txt | awk 'NR==2 {print $3}')

      if [[ $((NMSK * 5)) -ge $NLAB ]]; then
        VAL=$(cat /tmp/out.txt | awk 'NR==5 { print $4 - 1.0 }')
	MEAN_VAL=$(cat /tmp/out_mean.txt | awk 'NR==5 { print $2 - 1.0 }')
      else
        VAL="NA"
	MEAN_VAL="NA"
      fi

      echo $id,$LABSTR,$NLAB,$NMSK,$VAL, $MEAN_VAL | tee -a $OUTCSV

    done

  done
}

function analyze_density_hotspots()
{
  P=90
  OUTCSV=$ATLAS_DIR/tau_atlas/template_density_Tau_tangles_hotspots_${P}_tiles.csv
  TSEG=$ATLAS_DIR/map_hotspots_v2/hotspot_seg/hotspot_seg_combined.nii.gz
  TSEG_DIL=$ATLAS_DIR/map_hotspots_v2/hotspot_seg/hotspot_seg_combined_dilated.nii.gz

  TLAB=$ATLAS_DIR/histo_segs/hotspot_labels_itksnap.txt
  analyze_density $TSEG $TLAB $OUTCSV $P
}


function analyze_density_ba35()
{
  P=90
  OUTCSV=$ATLAS_DIR/tau_atlas//template_density_Tau_tangles_histoseg_unnormalized_${P}_tiles.csv
  TSEG=$ATLAS_DIR/histo_segs/template_hist_29_coarse/template_histo_seg_final.nii.gz
  TLAB=$ATLAS_DIR/histo_segs/template_histoseg_labels_itksnap.txt
  analyze_density $TSEG $TLAB $OUTCSV $P
}


main()
{

  SUBJ=$(cat $ATLAS_DIR/scripts/subj_tau_atlas.txt)
  N=$(cat $ATLAS_DIR/scripts/subj_tau_atlas.txt | wc -l)

  mkdir -p $WORK_DIR/dump
  #for((i=1;i<=$N;i++)); do

    #ID=$(cat $ATLAS_DIR/scripts/tau_atlas_subj.txt | head -n $i | tail -n 1 | awk '{print $1}')
    #vis_to_raw_space $ID
    #qsub -l h_vmem=15G,s_vmem=14.9G -cwd -o $WORK_DIR/dump -j y -N "tau_to_template_${ID}" $0 warp_tau_to_template $ID
    #warp_tau_to_template $ID
  #done

  #tau_map_average

  analyze_density_hotspots
  #analyze_density_ba35


}


main
