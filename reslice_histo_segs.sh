#/bin/bash
set -x -e

# atlas directory
ROOT=/project/hippogang_3/sravikumar/atlasPHG2019

source $ROOT/scripts/common.sh

#Template directory
TEMPLATE_ROOTDIR=$ROOT/atlasPHG_V3_2022
TDIR=$TEMPLATE_ROOTDIR/mst_multilabel
TEMPLATE_PHG=$TDIR/template_ncc/iter5/template_phg_thresh0.5.nii.gz
TEMPLATE_SRLM=$TDIR/template_ncc/iter5/template_srlm_thresh0.5.nii.gz

#Directory for histology work
HDIR=$ROOT/histo_segs

# Paul's directory with registrations for histology scans
TAU_ATLAS=/project/hippogang_2/pauly/tau_atlas

#Copy the histology annotations done by Sydney from a single folder into subject specific folder. Originally downloaded from box
function fn_copy_sydney_histo_annot()
{

 # Read input arguments
 id=${1?}
 block=${2?}

 id_PY=$(cat $ROOT/scripts/subj_histo_annot.txt | awk -v id=$id '{if ($1 == id) print $2}')

 WDIR=$ROOT/histo_segs/$id
 
 HISTO_SEG=$(find $ROOT/sydney_histo_annot_fromBox/ -name "${id_PY}_${block}_*wsrlm.nii.gz")
 if [[ -f $HISTO_SEG ]]; then

	mkdir -p $WDIR/histo_annot/${block}
        mv $HISTO_SEG $WDIR/histo_annot/${block}/${id}_${block}_histo_annot.nii.gz
 fi

}

# Copy histology warps from paul's folder
function fn_copy_paul_warps()
{
  mkdir -p $HDIR
  cp $TAU_ATLAS/manifest/blockface_src.txt $HDIR/
  N=$(cat $ROOT/scripts/subj_histo_annot.txt | wc -l)

  #Couldn't find the correctly saved old warps for 120126 in PY directory. Used the warps from the old approach that I had saved in my directory  
  for ((i=16;i<=$N;i++)); do

    id=$(cat $ROOT/scripts/subj_histo_annot.txt | head -n $i | tail -n 1)
    fn=$(echo $id | awk '{print $1}')
 
    fn_PY=$(cat $ROOT/scripts/subj_histo_annot.txt | awk -v id=$fn '{if ($1 == id) print $2}') 
    reg_mode=$(cat $ROOT/scripts/subj_histo_annot.txt | awk -v id=$fn '{if ($1 == id) print $3}')

    WDIR=$HDIR/$fn/histo_reg
    mkdir -p $WDIR

    read -r dummy blocks <<< "$(grep "^${fn_PY}" "$HDIR/blockface_src.txt")" 
    for block in $blocks; do

    echo $fn_PY $block
    mkdir -p $WDIR/$block

    #Check if this subject was registered using the old or new method developed by PY
    if [[ $reg_mode == 0 ]]; then
        echo "Warped using old registration method"
              
        # Old version of warp chain consists of many steps
        MAT=$(find $TAU_ATLAS/manual/${fn_PY}/hires_to_mold -name "${fn_PY}*_mri_hires_to_mold_affine.mat")
        #WARP=$(find $ROOT/histo_segs/rescue_hisseg_warp_chains/${fn_PY}/$block -name "${fn_PY}*_mri_invwarp_fx_hires_mv_mold.nii.gz")
        
        if [[ ! -d $TAU_ATLAS/brain_work/${fn_PY} ]];
        then 
          MOLD_TO_HIRES_ROOT_WARP=$(find $TAU_ATLAS/work/${fn_PY}/mri/ -name "${fn_PY}*_mri_rootwarp_fx_hires_mv_mold.nii.gz")
          #Per block affine
          AFFINE_FULL=$(find $TAU_ATLAS/work/${fn_PY}/bfreg/$block -name "${fn_PY}_${block}*_mri_to_bfvis_affine_full.mat")
        else
          MOLD_TO_HIRES_ROOT_WARP=$(find $TAU_ATLAS/brain_work/${fn_PY}/mri/ -name "${fn_PY}*_mri_rootwarp_fx_hires_mv_mold.nii.gz")
          AFFINE_FULL=$(find $TAU_ATLAS/brain_work/${fn_PY}/bfreg/$block -name "${fn_PY}_${block}*_mri_to_bfvis_affine_full.mat")
        fi

        #Output file
        HIRES_TO_MOLD_AFFINE=$WDIR/$block/${fn}_mri_hires_to_mold_affine.mat
        MOLD_TO_HIRES_WARP_ROOT=$WDIR/$block/${fn}_mri_rootwarp_fx_hires_mv_mold.nii.gz
        MRI_TO_BFVIS_AFFINE_FULL=$WDIR/$block/${fn}_${block}_mri_to_bfvis_affine_full.mat

        cp $MAT $HIRES_TO_MOLD_AFFINE
        cp $MOLD_TO_HIRES_ROOT_WARP $MOLD_TO_HIRES_WARP_ROOT
        cp $AFFINE_FULL $MRI_TO_BFVIS_AFFINE_FULL

        WARP_CHAIN_BFVIS_TO_HIRES_RAW="$MOLD_TO_HIRES_WARP_ROOT,64 $HIRES_TO_MOLD_AFFINE,-1 $MRI_TO_BFVIS_AFFINE_FULL,-1"
        echo $WARP_CHAIN_BFVIS_TO_HIRES_RAW > $WDIR/$block/chain_bfvis_to_hires_${block}.txt
        
        WARP_CHAIN_HIRES_TO_BFVIS="$MRI_TO_BFVIS_AFFINE_FULL $HIRES_TO_MOLD_AFFINE $MOLD_TO_HIRES_WARP_ROOT,-64"
        echo $WARP_CHAIN_HIRES_TO_BFVIS > $WDIR/$block/chain_hires_to_bfvis_${block}.txt	

      elif [[ $reg_mode == 1 ]]; then

        echo "Using warps from new approach"
        HIRES_TO_BFVIS_WARP_FULL=$(find $TAU_ATLAS/work/$fn_PY/bfreg/$block -name "${fn_PY}_${block}_hires_mri_to_bfvis_full_warp.nii.gz")
        HIRES_TO_BFVIS_INVWARP_FULL=$(find $TAU_ATLAS/work/$fn_PY/bfreg/$block -name "${fn_PY}_${block}_hires_mri_to_bfvis_full_invwarp.nii.gz")

        cp $HIRES_TO_BFVIS_WARP_FULL $WDIR/$block/${fn}_${block}_hires_mri_to_bfvis_full_warp.nii.gz
        cp $HIRES_TO_BFVIS_INVWARP_FULL $WDIR/$block/${fn}_${block}_hires_mri_to_bfvis_full_invwarp.nii.gz

        WARP_CHAIN_BFVIS_TO_HIRES_RAW="$WDIR/$block/${fn}_${block}_hires_mri_to_bfvis_full_invwarp.nii.gz"
        echo $WARP_CHAIN_BFVIS_TO_HIRES_RAW > $WDIR/$block/chain_bfvis_to_hires_${block}.txt

        WARP_CHAIN_HIRES_TO_BFVIS="$WDIR/$block/${fn}_${block}_hires_mri_to_bfvis_full_warp.nii.gz"
        echo $WARP_CHAIN_HIRES_TO_BFVIS > $WDIR/$block/chain_hires_to_bfvis_${block}.txt
      fi
      
	   
    done

  done
}

# To register ex vivo hotspots to histology space. Use the ex vivo hotspot image generated for in vivo analysis
# For paper figure showing subfields corresponding to hotspot in subject space
function fn_map_hotspots_to_histo()
{

   export LD_LIBRARY_PATH=/project/hippogang_3/sravikumar/atlasPHG2019/pkgs/libc/:/project/hippogang_3/sravikumar/packages/glibc-2.14/lib:/project/hippogang_3/sravikumar/atlasPHG2019/pkgs:/opt/python/lib

 # Read input arguments
 id=${1?}
 block=${2?}

 ID_PY=$(cat $ROOT/scripts/subj_map_hotspots.txt | awk -v id=$id '{if ($1 == id) print $2}')
 WDIR=$HDIR/$id
 WORKDIR=$WDIR/hotspots

 # Directory to store warped hotspots
 mkdir -p $WORKDIR

 HIRES_MRI_MANUAL_PHGSEG_AFFINE=$(find $ROOT/inputs/ -name "${id}*_raw_to_axisalign.mat")
 # Raw to histo space transformations - from Paul 
 BFVIS_HIRES_MRI_RESIDUAL_T0_BF_WARP=$WDIR/histo_reg/${id}_${block}_mri_residual_to_bfvis_warp.nii.gz
 MRI_TO_BFVIS_INIT_RIGID=$WDIR/histo_reg/${id}_${block}_mri_to_bfvis_init_rigid.mat
 BFVIS_HIRES_MRI_RESIDUAL_TO_BF_AFFINE=$WDIRhisto_reg//${id}_${block}_mri_residual_to_bfvis_affine.mat
 MRI_TO_BFVIS_AFFINE_FULL=$WDIR/histo_reg/${id}_${block}_mri_to_bfvis_affine_full.mat
 HIRES_TO_MOLD_AFFINE=$WDIR/histo_reg/${id}_mri_hires_to_mold_affine.mat
 MOLD_TO_HIRES_INV_WARP=$WDIR/histo_reg/${id}_mri_invwarp_fx_hires_mv_mold.nii.gz
  
if [[ -f $MRI_TO_BFVIS_INIT_RIGID ]];then
	HIRES_TO_BFVIS_WARP_FULL="$BFVIS_HIRES_MRI_RESIDUAL_T0_BF_WARP $BFVIS_HIRES_MRI_RESIDUAL_TO_BF_AFFINE $MRI_TO_BFVIS_INIT_RIGID $HIRES_TO_MOLD_AFFINE $MOLD_TO_HIRES_INV_WARP"
else
	HIRES_TO_BFVIS_WARP_FULL="$MRI_TO_BFVIS_AFFINE_FULL $HIRES_TO_MOLD_AFFINE $MOLD_TO_HIRES_INV_WARP"
fi


# Template to ex vivo subject space
SDIR=$TDIR/final/${id}
WARP_CHAIN_TEMPLATE_TO_UNWARP=$(cat "$SDIR/chain_template_to_unwarped.txt")
WARP_CHAIN_TEMPLATE_TO_RAW="$ROOT/preproc/$id/unwarp/${id}_unwarp.mat,-1 $WARP_CHAIN_TEMPLATE_TO_UNWARP"

 #Reference_histology_space
 #HISTO_ANNOT_SPLAT_MANIFEST=$TAU_ATLAS/work/$ID_PY/historeg/$block/splat/${ID_PY}_${block}_annot_splat_manifest.txt
 #HISTO_ANNOT_SPLAT_PATTERN=$TAU_ATLAS/work/$ID_PY/historeg/$block/splat/${ID_PY}_${block}_nissl_mrilike_splat_%s.nii.gz

 mkdir -p $WDIR/refspace
 REF_HISTO_SPACE=$WDIR/refspace/${id}_${block}_histo_refspace.nii.gz
 BFVIS_MRILIKE=$TAU_ATLAS/work/$ID_PY/bfreg/${block}/${ID_PY}_${block}_bfvis_mrilike.nii.gz
 cp $BFVIS_MRILIKE $REF_HISTO_SPACE

<<"OLD"
 if [[ ! -f $REF_HISTO_SPACE ]]; then

        STAGE=voliter-20
        #if [[ -f $HISTO_ANNOT_SPLAT_MANIFEST ]]; then
          local OUT="$(printf $HISTO_ANNOT_SPLAT_PATTERN $STAGE)"
        #fi

         cp $OUT $REF_HISTO_SPACE
 fi
OLD

 #Hotspot segmentation in ex vivo atlas space
 HOTSPOT_SEG="$ROOT/map_hotspots_v2/hotspot_seg/hotspot_seg_combined.nii.gz"

 #Output location
 HISTO_NISSL_SPLAT_HOTSPOT_SEG=$WORKDIR/${id}_${block}_nissl_splat_tau_hotspot.nii.gz
 
 #hires mri warped to bfvis space
 cp $TAU_ATLAS/work/$ID_PY/bfreg/${block}/${ID_PY}_${block}_hires_mri_to_bfvis_warped.nii.gz $WORKDIR/${id}_${block}_hires_mri_to_bfvis_warped.nii.gz

 # This is a special case. the raw scan used in the histology recon has the wrong header info.
 # The corrected raw scan was used during atlas construction. Need to inlclude an additional regitration step
 # between these two raw scans. A smilar situation is likely for HNL11-15 when we come to it.
 if [[ "${id}" == "INDD116748-R" ]];then


   HISTO_RAW=$(find $ROOT/pmatlas_raw/ -name "*116748*_wrongheader.nii.gz")
   ATLAS_RAW=$(find $ROOT/pmatlas_raw/ -name "*116748*_stitch_20160314.nii")
   MAT=$WDIR/histo_reg/${id}_mri_atlas_to_histo_raw.mat

<<"DONE"
   if [[ ! -f $MAT ]]; then
        $GREEDY_HOME/greedy -d 3 -a \
        -i $HISTO_RAW $ATLAS_RAW \
        -dof 6 -n 100x50x0 -m NCC 4x4x4 \
        -o $MAT
   fi
DONE
   WARP_CHAIN_RAW_TO_HISTO="$HIRES_TO_BFVIS_WARP_FULL $MAT"

 else

  WARP_CHAIN_RAW_TO_HISTO="$HIRES_TO_BFVIS_WARP_FULL"

 fi

 WARP_CHAIN_TEMPLATE_TO_HISTO="$WARP_CHAIN_RAW_TO_HISTO $WARP_CHAIN_TEMPLATE_TO_RAW"

   $GREEDY_HOME/greedy -d 3 \
  -rf $REF_HISTO_SPACE -ri LABEL 0.04mm \
  -rm $HOTSPOT_SEG $HISTO_NISSL_SPLAT_HOTSPOT_SEG \
  -r $WARP_CHAIN_TEMPLATE_TO_HISTO

}

function fn_map_histo_segs_to_hiresmri() 
{

 # Read input arguments
 id=${1?}
 block=${2?}
 
 seg_mode=${3?}
 
 id_PY=$(cat $ROOT/scripts/subj_histo_annot.txt | awk -v id=$id '{if ($1 == id) print $2}')
 
 WDIR=$ROOT/histo_segs/$id
 WORK=$WDIR/wholemri_${seg_mode}/${block}
 mkdir -p $WORK

 # Raw to axisalign space
 HIRES_MRI_MANUAL_PHGSEG_AFFINE=$(find $ROOT/inputs/ -name "${id}*_raw_to_axisalign.mat")

 if [[ ! -f $HIRES_MRI_MANUAL_PHGSEG_AFFINE ]]; then
        HIRES_MRI_MANUAL_PHGSEG_AFFINE=$ROOT/scripts/manual/raw_to_axisalign_transforms/${id}_raw_to_axisalign.mat
 fi

 WARP_CHAIN_BFVIS_TO_HIRES=$(cat $WDIR/histo_reg/$block/chain_bfvis_to_hires_${block}.txt)
 
 # Reference axisalign space
 MRI_AXISALIGN=$(find $ROOT/inputs/ -name "${id}*_axisalign_img.nii.gz")

 #Use the same raw scans used in histo recon. Minor differences in scans lead to errors 
 MRI_RAW=$WDIR/${id}_mri_hires_historecon.nii.gz
 if [[ ! -f $MRI_RAW ]]; then
   MRI_RAW_PAUL=$TAU_ATLAS/input/${id_PY}/hires_mri/${id_PY}_mri_hires.nii.gz
   cp -u $MRI_RAW_PAUL $MRI_RAW
 fi

 #Histology annotation done by Sydney
 HISTO_ANNOT=$WDIR/histo_annot/${block}/${id}_${block}_histo_annot.nii.gz
 
 #Intermediate/temporary files
 SLICES=$WORK/slices.nii.gz
 SLICESBG=$WORK/slices_bg.nii.gz
 INSIDE_REGION=$WORK/mask_inside_region.nii.gz

if [[ $seg_mode == "coarse" ]];then

	#Combine ERC sufields, BA35 subfields, BA36 etc
	c3d $HISTO_ANNOT -as SEG \
	-retain-labels 4 5 13 14 15 20 21 22 23 26 37 42 44 45 46 51 53 56 57 59 61 62 63 \
	-replace  44 42 56 53 57 53 62 61 63 61 -o $WORK/${id}_${block}_histo_annot_subseg1.nii.gz \
	-push SEG -thresh 6 9 6 0 -as BA35 -o $WORK/${id}_${block}_histo_annot_subsegba35.nii.gz \
	-push SEG -retain-labels 11 12 17 -thresh 11 17 11 0 -as BA36 \
	-o $WORK/${id}_${block}_histo_annot_subsegba36.nii.gz \
	-push SEG -thresh 27 36 27 0 -as ERC -o $WORK/${id}_${block}_histo_annot_subsegerc.nii.gz 

  HISTO_ANNOT=$WORK/${id}_${block}_histo_annot_coarseseg.nii.gz

	c3d $WORK/${id}_${block}_histo_annot_subseg*.nii.gz \
	-accum -add -endaccum -o $HISTO_ANNOT

fi

 #if [[ ! -f $WORK/slices_bg_to_wholemri.nii.gz ]]; then

	    c3d $HISTO_ANNOT -as X -thresh 1 inf 1 0 -dilate 1 200x200x0 -push X -times -o $SLICES
	    c3d $HISTO_ANNOT -as X -thresh 1 inf 1 0 -dilate 1 200x200x0 -push X -thresh 1 inf 1 0 -times \
	    -as Y -dilate 1 600x600x0 -push Y -scale -1 -add \
	    -replace 1 70 $SLICES -add \
	    -trim -0vox -o $SLICESBG
  #fi

	  # Create a mask between the first and last slice - what should be interpolated
	  c3d $HISTO_ANNOT -thresh 1 inf 1 0 -dilate 1 10x10x0 -trim 0mm -thresh -inf inf 1 0 \
	    -pad 2x2x2 2x2x2 0 -o $INSIDE_REGION

 if [[ "${id_PY}" == "INDD116748" ]];then

	  HISTO_RAW=$MRI_RAW #$(find $ROOT/pmatlas_raw/ -name "*116748*_wrongheader*.nii.gz")
   	#ATLAS_RAW=$(find $ROOT/pmatlas_raw/ -name "*116748*_stitch_20160314.nii")
   	ATLAS_RAW=$ROOT/preproc/${id}/${id}_raw_n4clip.nii.gz
	  ATLAS_TO_HISTO_RAW_MAT=$WDIR/histo_reg/${id}_mri_atlas_to_histo_raw.mat

   	if [[ ! -f $ATLAS_TO_HISTO_RAW_MAT ]]; then
   	     $GREEDY_HOME/greedy -d 3 -a \
   	     -i $HISTO_RAW $ATLAS_RAW \
   	     -dof 6 -n 100x50x0 -m NCC 4x4x4 \
   	    -o $ATLAS_TO_HISTO_RAW_MAT
   	fi

	  ln -sf $ATLAS_RAW $WDIR/${id}_mri_hires_atlas.nii.gz

    $GREEDY_HOME/greedy -d 3 \
    -rf $ATLAS_RAW \
    -ri LABEL 0.2vox \
    -rm $SLICESBG $WORK/slices_bg_to_wholemri.nii.gz \
    -rm $INSIDE_REGION $WORK/inside_mask_to_wholemri.nii.gz \
    -rm $SLICES $WORK/slices_to_wholemri.nii.gz \
    -r $ATLAS_TO_HISTO_RAW_MAT,-1 $WARP_CHAIN_BFVIS_TO_HIRES

 else	
	   $GREEDY_HOME/greedy -d 3 \
	   -rf $MRI_RAW \
	   -ri LABEL 0.2vox \
	   -rm $SLICESBG $WORK/slices_bg_to_wholemri.nii.gz \
	   -rm $INSIDE_REGION $WORK/inside_mask_to_wholemri.nii.gz \
	   -rm $SLICES $WORK/slices_to_wholemri.nii.gz \
	    -r $WARP_CHAIN_BFVIS_TO_HIRES

fi
 
## Adding this step since I was having issues with INDD119294. Once I map slices_bg to wholemri space,
# The background has a value of 5 (not zero) which causes errors in the step combining the segmentations
c3d $WORK/inside_mask_to_wholemri.nii.gz $WORK/slices_bg_to_wholemri.nii.gz -times -o $WORK/slices_bg_to_wholemri_masked.nii.gz

}

function combine_blocks()
{

  ID=${1?}
  SPACE=wholemri

  seg_mode=${2?}

  # Combine the slices from all the blocks
  WDIR=$ROOT/histo_segs/$ID
  WORK=$WDIR/wholemri_${seg_mode}/

  # Link to whole-brain MRI
  MRI_RAW=$WDIR/${ID}_mri_hires_historecon.nii.gz  
 
  # Combined slices
  COMB_SLICES_NOBG=$WORK/${ID}_combine_slices_to_${SPACE}.nii.gz
  COMB_SLICES=$WORK/${ID}_combine_slices_bg_to_${SPACE}.nii.gz
  RESULT=$WORK/${ID}_histoseg_interp_${SPACE}.nii.gz
  MASK=$WORK/${ID}_block_mask_${SPACE}.nii.gz
  RESULT_MASKED=$WORK/${ID}_histoseg_interp_masked_${SPACE}.nii.gz


<<"SKIP1"
  # Add segmentations
  c3d $WORK/*/slices_to_${SPACE}.nii.gz \
    -foreach -thresh 1 inf 1 0 -endfor -accum -add -endaccum -thresh 1 1 1 0 \
    -popas MASK \
    $WORK/*/slices_to_${SPACE}.nii.gz \
    -accum -add -endaccum -push MASK -times \
    -o $COMB_SLICES_NOBG
SKIP1

  # Add segmentations
  c3d $WORK/*/slices_bg_to_${SPACE}_masked.nii.gz \
    -foreach -thresh 1 inf 1 0 -endfor -accum -add -endaccum -thresh 1 1 1 0 \
    -popas MASK \
    $WORK/*/slices_bg_to_${SPACE}_masked.nii.gz \
    -accum -add -endaccum -push MASK -times \
    -o $COMB_SLICES

  # Interpolate using distance transform
  c3d -verbose $COMB_SLICES \
    -replace 0 100 -split \
    -foreach -sdt -clip 0 inf -smooth-fast 0.2mm -scale -1 -endfor \
    -scale 0 -shift -1e10 -merge \
    -replace 70 0 -o $RESULT

  # Combine block binary masks to create a mask of covered regions - still need this step after interpolating. 
  c3d $WORK/*/inside_mask_to_wholemri.nii.gz \
    -accum -add -endaccum -thresh 1 1 1 0 \
    -o $MASK \
    $RESULT -times -o $RESULT_MASKED

}

function warp_subject_to_template()
{

   ID=${1?}

   #Reference image for new cases is the existing template
   REF_TEMPLATE_CSMED=$TDIR/template_ncc/iter5/template_csmed.nii.gz
   REF_TEMPLATE_CSLAT=$TDIR/template_ncc/iter5/template_cslat.nii.gz
   REF_TEMPLATE_SRLM=$TDIR/template_ncc/iter5/template_srlm.nii.gz
   REF_TEMPLATE_IMG=$TDIR/template_ncc/iter5/template_img.nii.gz
   REF_TEMPLATE_PHG=$TDIR/template_ncc/iter5/template_phg_thresh0.5.nii.gz

   #Input images
   CSMEDMOV=$ROOT/preproc/$ID/unwarp/multi_label/${ID}_unwarp_csmedsegshape.nii.gz
   CSLATMOV=$ROOT/preproc/$ID/unwarp/multi_label/${ID}_unwarp_cslatsegshape.nii.gz
   SRLMMOV=$ROOT/preproc/$ID/unwarp/${ID}_unwarp_srlmseg_sr.nii.gz
   PHGMOV=$ROOT/preproc/$ID/unwarp/${ID}_unwarp_phgsegshape_singlelabel.nii.gz
 
   ODIR=$HDIR/${ID}/reg_to_templateAD
   mkdir -p $ODIR
   TAG="${ID}_to_template"

   MAT_MOMENTS=$ODIR/${TAG}_moments.mat
   AFFINE=$ODIR/${TAG}_affine.mat
   WARPROOT=$ODIR/${TAG}_rootwarp.nii.gz
   WARP=$ODIR/${TAG}_warp.nii.gz

   $GREEDY_HOME/greedy -d 3 \
   -i $REF_TEMPLATE_PHG $PHGMOV \
   -moments -m NCC 4x4x4 \
   -o $MAT_MOMENTS
  
   $GREEDY_HOME/greedy -d 3 \
   -w 1000 -i $REF_TEMPLATE_CSMED $CSMEDMOV \
   -w 2000 -i $REF_TEMPLATE_CSLAT $CSLATMOV \
   -w 5000 -i $REF_TEMPLATE_SRLM $SRLMMOV \
   -n 100x100 -m SSD \
   -ia $MAT_MOMENTS -a -o $AFFINE

   # Perform registration - use intensity along with boundary strength - took out -w 0.00006 -i $IMGREF $IMGMOV
   # Increase number of iterations from 100x50x40x20 - increase smoothing from 0.6/0.1
    $GREEDY_HOME/greedy -d 3 \
    -w 1000 -i $REF_TEMPLATE_CSMED $CSMEDMOV \
    -w 2000 -i $REF_TEMPLATE_CSLAT $CSLATMOV \
    -w 5000 -i $REF_TEMPLATE_SRLM $SRLMMOV \
    -n 100x100x50x50 \
    -s 0.6mm 0.1mm -e 0.5 -float \
    -it $AFFINE \
    -sv -wp 0 -oroot $WARPROOT \
    -o $WARP
#was 0.8mm and 50 at 2nd level

}

function process_seg()
{
  # This is the subject ID
  id=${1?}

  # Registration mode for which this is being generated
  mode=${2?}

  # Coarse or fine segmentation
  seg_mode=${3?}

  # Get the histology segmentation filename
  HIST_SEG=$HDIR/$id/wholemri_${seg_mode}/${id}_histoseg_interp_masked_wholemri.nii.gz
  
  #Work directory
  WORK=$TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}

  # Does Sydney's histology segmentation exist?
  if [[ -f $HIST_SEG ]]; then
    SEG=$HIST_SEG
  else
    echo "MISSING SEGMENTATION!"
    exit -1
  fi

  #The 2022 atlas included in dissertation doesn't include many of the subjects with histology annotations. Map these segmentations directly to the template. Atlas built in Aug 2022 (V3) includes these cases to avoid this
  # instead of using atlas transforms
  if [[ -z $(grep ${id} $ROOT/scripts/subj_atlas_2022_v2.txt) ]]; then
	  echo "Subject not included in atlas"
	  ID_SIDE=$id
    warp_subject_to_template ${ID_SIDE}
	  ATLAS=0
  else 
    # The subjects ids used during atlas construction include the side
    ID_SIDE=$id
    ATLAS=1
  fi
 
  # The histology segmentations are in raw image space. Used the raw scan used by Paul. Should be the 
  # same as the raw scan used by Sadhana in the atlas but might be slight differences in header info. 
  #This is def the case for INDD116748
  IMG=$ROOT/preproc/$ID_SIDE/${ID_SIDE}_raw_n4clip.nii.gz

  # Check that we have all files
  if [[ ! -f $IMG ]]; then
    echo "Missing raw image"
    exit
  fi

  # For warping the histology segmentation to axisalign space - this is where manual segmentations are done
  # Raw to axisalign space
  HIRES_MRI_MANUAL_PHGSEG_AFFINE="$ROOT/preproc/${ID_SIDE}/${ID_SIDE}_transform_to_axisalign.mat"

  # Reference axislaign and unwarped space
  #MRI_AXISALIGN="$ROOT/inputs/${ID_SIDE}/${ID_SIDE}_axisalign_img.nii.gz"
  MRI_AXISALIGN=$ROOT/preproc/${ID_SIDE}/${ID_SIDE}_axisalign_img.nii.gz
	 
  #Distance map boundary segmentations in axisalign space. Needed to compaire volumetric vs unfolded registration for 
  # MICCAI 2021 submission
  #UNFOLD_GM=$ROOT/preproc/$ID_SIDE/${ID_SIDE}_unfolded_gm.nii.gz
  #$C3D_HOME/c3d $ROOT/inputs/${ID_SIDE}/${ID_SIDE}_labelmap.nii.gz -retain-labels 1 -o $UNFOLD_GM

<<"DONTNEED" 
  # There often seems to be a difference between the raw space used for the atlas construction and histology reconstruction
  # Include a registration between the two raw scans
   HISTO_RAW=$HDIR/${id}/${id}_mri_hires_historecon.nii.gz
   ATLAS_RAW=$IMG

   MAT=$WDIR/histo_reg/${id}_mri_atlas_to_histo_raw.mat

   if [[ ! -f $MAT ]]; then
        greedy -d 3 -a \
        -i $HISTO_RAW $ATLAS_RAW \
        -dof 6 -n 100x50x0 -m NCC 4x4x4 \
        -o $MAT
   fi

  # Make sure we have the correct header
  #c3d $IMG $SEG -copy-transform -o $TMPDIR/${id}_seg.nii.gz
DONTNEED

  # We also need to create a mask of the segmentation
  $C3D_HOME/c3d $SEG -thresh 1 inf 1 0 -o $TMPDIR/${id}_mask.nii.gz

  if [[ $ATLAS == 1 ]]; then
  # Set up the warpchains and working directories for the different
  # modes that are available
  local WDIR WARPCHAIN
  case "$mode" in
    ncc)
      WDIR=$WORK
      sed -i 's#/data/jux/#/project/hippogang_3/#gi' $TDIR/final/${ID_SIDE}/chain_raw_to_template.txt
      WARPCHAIN=$(cat $TDIR/final/${ID_SIDE}/chain_raw_to_template.txt)
      ;;
    aff_ncc)
      WDIR=$WORK/$mode
      WARPCHAIN="\
        $TEMPLATE_AFF_DIR/warp_${ID_SIDE}.nii.gz \
        $TDIR/template_aff_ncc/input/$ID_SIDE/proc_affine_${ID_SIDE}.mat \
	$ROOT/preproc/$ID_SIDE/unwarp/${ID_SIDE}_unwarp.mat"
      ;;
    mst)
      WDIR=$WORK/$mode
      WARPCHAIN="$(cat $TDIR/paths/$ID_SIDE/final/chain_unwarp_to_final.txt) \
	$ROOT/preproc/$ID_SIDE/unwarp/${ID_SIDE}_unwarp.mat"
      ;;
    gshoot)
      WDIR=$WORK/$mode
      WARPCHAIN="\
        $TDIR/gshoot/${ID_SIDE}/iter_4/shooting_warp.nii.gz \
        $TDIR/gshoot/${ID_SIDE}/iter_4/target_to_root_procrustes.mat,-1 \
	$ROOT/preproc/$ID_SIDE/unwarp/${ID_SIDE}_unwarp.mat"
      ;;
  esac

 elif [[ $ATLAS == 0 ]] ; then
	WDIR=$WORK
	WARPCHAIN="$HDIR/${ID_SIDE}/reg_to_templateAD/${ID_SIDE}_to_template_warp.nii.gz $HDIR/${ID_SIDE}/reg_to_templateAD/${ID_SIDE}_to_template_affine.mat $ROOT/preproc/$ID_SIDE/unwarp/${ID_SIDE}_unwarp.mat"
 fi

<<"DONE"
  # Create the working directory
  WORK=$WDIR/histoseg_mapped_totemplate
  mkdir -p $WORK
  
  TEMPLATE_SEG=$TDIR/template_ncc/iter5/template_phg_thresh0.5.nii.gz

  # Transform the segmentation into the template space
  $GREEDY_HOME/greedy -d 3 \
    -rf $TDIR/template_ncc/iter5/template_img.nii.gz \
    -rm $IMG $WORK/${id}_wholemri_to_template.nii.gz \
    -ri LABEL 0.24vox \
    -rm $TMPDIR/${id}_mask.nii.gz $WORK/${id}_histseg_to_template_mask.nii.gz \
    -rm $SEG $WORK/${id}_histseg_to_template.nii.gz \
    -r $WARPCHAIN 

 c3d $WORK/${id}_histseg_to_template_mask.nii.gz $TEMPLATE_SEG -times \
 -o  $WORK/${id}_histseg_to_template_mask_phgmasked.nii.gz 

  c3d $WORK/${id}_histseg_to_template.nii.gz $TEMPLATE_SEG -times \
 -o  $WORK/${id}_histseg_to_template_phgmasked.nii.gz
DONE

 # Warp the histology segmentation to axisalign space
 ADIR=$HDIR/$id/axisalign/
 mkdir -p $ADIR

 $GREEDY_HOME/greedy -d 3 \
 -rf $MRI_AXISALIGN \
 -ri LABEL 0.24vox \
 -rm $SEG $ADIR/${id}_histoseg_interp_masked_wholemri_to_axisalign.nii.gz \
 -r $HIRES_MRI_MANUAL_PHGSEG_AFFINE

 AXISALIGN_PHGSEG=$ROOT/preproc/${id}/${id}_axisalign_phgsegshape_singlelabel.nii.gz 
 $C3D_HOME/c3d $ADIR/${id}_histoseg_interp_masked_wholemri_to_axisalign.nii.gz $AXISALIGN_PHGSEG -times \
 -o $ADIR/${id}_histoseg_interp_masked_wholemri_to_axisalign_corrected.nii.gz 

<<"SKIP"
# Warp gray matter label map from axislaign space to template - for MICCAI unfolding paper
 $GREEDY_HOME/greedy -d 3 \
 -rf $TDIR/template_ncc/iter5/template_img.nii.gz \
 -ri LABEL 0.24vox \
 -rm $UNFOLD_GM $WDIR/unfolded_histo_seg/${id}_histseg_unfoldmask_to_template.nii.gz \
 -r $WARPCHAIN $HIRES_MRI_MANUAL_PHGSEG_AFFINE,-1
SKIP

}

function combine_histo_subj()
{

  mode=${1?}

  seg_mode=${2?}

  WORK=$TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}
  local WDIR
    
  # Work directory
  if [[ $mode == "ncc" ]]; then
    WDIR=$WORK
  else
    WDIR=$WORK/$mode
  fi

   mkdir -p $WORK

 #Don't include labels 1 and 2 (coarse seg labels) in the sum. 
 c3d $WDIR/histoseg_mapped_totemplate/*histseg_to_template_phgmasked.nii.gz -foreach -thresh 4 inf 1 0 -endfor -accum -add -endaccum -o $WDIR/sum.nii.gz

 mkdir -p $WDIR/hitsum
 mkdir -p $WDIR/probmap

 #added 13 (Jan 3 2023)
 if [[ $seg_mode == "fine" ]]; then

	  for i in 5 6 7 8 9 11 12 13 17 20 21 22 23 26 27 30 31 32 33 34 35 36 37 42 44 45 46 51 53 56 57 59 61 62 63; do
	    c3d $WDIR/histoseg_mapped_totemplate/*histseg_to_template_phgmasked.nii.gz \
	      -foreach -thresh $i $i 1 0 -endfor \
	      -accum -add -endaccum \
	      -o $WDIR/hitsum/hitsum_${i}.nii.gz \
	      $WDIR/sum.nii.gz -reciprocal -times \
	      -o $WDIR/probmap/probmap_${i}.nii.gz
	  done
  elif [[ $seg_mode == "coarse" ]]; then

    for i in 5 6 11 13 20 21 22 23 26 27 37 42 45 46 51 53 59 61; do
      c3d $WDIR/histoseg_mapped_totemplate/*histseg_to_template_phgmasked.nii.gz \
        -foreach -thresh $i $i 1 0 -endfor \
        -accum -add -endaccum \
        -o $WDIR/hitsum/hitsum_${i}.nii.gz \
        $WDIR/sum.nii.gz -reciprocal -times \
        -o $WDIR/probmap/probmap_${i}.nii.gz
      done

  fi

  # Vote - no masking for now
  c3d $WDIR/sum.nii.gz -thresh 0 inf 0.1 0 \
    $WDIR/probmap/probmap_*.nii.gz -vote \
    -o $WDIR/vote_raw.nii.gz
   
  c3d $TEMPLATE_PHG $TEMPLATE_SRLM -add -o $WDIR/vote_hfdb.nii.gz
  
   c3d $WDIR/vote_hfdb.nii.gz -as A -thresh 1 1 1 0 $WDIR/vote_raw.nii.gz -times -push A -thresh 2 2 5 0 -add -o $WDIR/vote_masked.nii.gz

<<"SKIP"
  # Generate an entropy map and a masked entropy map
  c3d $WDIR/probmap/probmap_*.nii.gz \
    -foreach -dup -shift 0 -log -scale 1.442695 -clip -1000 1000 -times -endfor \
    -accum -add -endaccum -scale -1 -replace nan 0 \
    -o $WDIR/entropy.nii.gz

  c3d $WDIR/entropy.nii.gz $TEMPLATE_PHG \
    $TEMPLATE_PHG -thresh 1 1 0 -1 -add \
    -o $WDIR/entropy_masked.nii.gz
SKIP

}

function all_pair_histo_dice()
{

  seg_mode=${1?}
  local WDIR

  WORK=$TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}

  # Work directory
  WDIR=$WORK/histoseg_mapped_totemplate

  # List the available segmentations
  local IDS=($(for i in $(ls $WDIR/*_histseg_to_template.nii.gz); do basename $i _histseg_to_template.nii.gz; done))

  # Create a Dice directory
  mkdir -p $WDIR/overlap
  rm -rf $WDIR/gdsc.txt

  for ((i=0;i<${#IDS[*]};i++)); do
    for ((j=$((i+1));j<${#IDS[*]};j++)); do

      local OVLFILE=$WDIR/overlap/overlap_${IDS[i]}_${IDS[j]}.txt

      c3d \
        $TEMPLATE_PHG -popas M \
        $WDIR/${IDS[i]}_histseg_to_template.nii.gz -push M -times -as A \
        $WDIR/${IDS[j]}_histseg_to_template.nii.gz -push M -times -as B \
        -foreach -thresh 1 inf 1 0 -endfor -times -popas X \
        -push A -push X -times \
        -push B -push X -times \
        -label-overlap > $OVLFILE

      echo "$mode ${IDS[i]} ${IDS[j]} $(cat $OVLFILE | awk 'NR==3 {print $3}')" >> $WDIR/gdsc.txt

    done
  done

<<"SKIP"
  # Create a Dice directory for the segmentations masked over the unfolded regions only
  mkdir -p $WORK/unfolded_histo_seg/overlap_unfold
  UNFOLD_DIR=$WORK/unfolded_histo_seg
  rm -rf $UNFOLD_DIR/overlap/gdsc_unfold.txt

  for ((i=0;i<${#IDS[*]};i++)); do
    for ((j=$((i+1));j<${#IDS[*]};j++)); do

      local OVLFILE=$WORK/unfolded_histo_seg/overlap_unfold/overlap_${IDS[i]}_${IDS[j]}.txt

      $C3D_HOME/c3d \
        $TEMPLATE_PHG -popas M \
        $WDIR/${IDS[i]}_histseg_to_template.nii.gz $UNFOLD_DIR/${IDS[i]}_histseg_unfoldmask_to_template.nii.gz -times -push M -times -as A \
        $WDIR/${IDS[j]}_histseg_to_template.nii.gz $UNFOLD_DIR/${IDS[j]}_histseg_unfoldmask_to_template.nii.gz -times -push M -times -as B \
        -foreach -thresh 1 inf 1 0 -endfor -times -popas X \
        -push A -push X -times \
        -push B -push X -times \
        -label-overlap > $OVLFILE

      echo "$mode ${IDS[i]} ${IDS[j]} $(cat $OVLFILE | awk 'NR==3 {print $3}')" >> $UNFOLD_DIR/overlap/gdsc_unfold.txt
       # 6 11 27 42 46
      echo "$mode ${IDS[i]} ${IDS[j]} $(cat $OVLFILE | awk '$1==6  {print $1}') $(cat $OVLFILE | awk '$1==6 {print $4}')" >> $UNFOLD_DIR/overlap/dsc_6_unfold.txt
      echo "$mode ${IDS[i]} ${IDS[j]} $(cat $OVLFILE | awk '$1==11 {print $1}') $(cat $OVLFILE | awk '$1==11 {print $4}')" >> $UNFOLD_DIR/overlap/dsc_11_unfold.txt
      echo "$mode ${IDS[i]} ${IDS[j]} $(cat $OVLFILE | awk '$1==27 {print $1}') $(cat $OVLFILE | awk '$1==27 {print $4}')" >> $UNFOLD_DIR/overlap/dsc_27_unfold.txt
      echo "$mode ${IDS[i]} ${IDS[j]} $(cat $OVLFILE | awk '$1==42 {print $1}') $(cat $OVLFILE | awk '$1==42 {print $4}')" >> $UNFOLD_DIR/overlap/dsc_42_unfold.txt
      echo "$mode ${IDS[i]} ${IDS[j]} $(cat $OVLFILE | awk '$1==46 {print $1}') $(cat $OVLFILE | awk '$1==46 {print $4}')" >> $UNFOLD_DIR/overlap/dsc_46_unfold.txt

      done
  done
SKIP

}

function mrf()
{

  # The regularization weight
  WEIGHT=${1?}

  seg_mode=${2?}

  # The weight converted to something that can be in a directory name
  DIRNAME=$(echo $WEIGHT | awk '{printf "%06d\n", $1 * 100000}')
  WORK=$TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}
  WDIR=$WORK/mrf_vote/mrf_vote_$DIRNAME

  # Scaling
  SCALE=$(ls $WORK/histoseg_mapped_totemplate/*_histseg_to_template_phgmasked.nii.gz | wc -l | awk '{print 1.0 / $1}')

  # Create the directory
  mkdir -p $WDIR

  # Perform the graph cut optimization
  $C3D_HOME/c3d $WORK/probmap/probmap_*.nii.gz \
  $TEMPLATE_PHG -popas M \
  -foreach -push M -thresh 0 0 NaN 1 -times -endfor \
  -verbose -vote-mrf VA $WEIGHT \
  -o $WDIR/mrf_seg_masked.nii.gz

<<"SKIP"
  c3d \
    $WORK/hitsum/hitsum_*.nii.gz \
    $TEMPLATE_PHG -popas M \
    -foreach -scale $SCALE -push M -thresh 0 0 NaN 1 -times -endfor \
    -verbose -vote-mrf VA $WEIGHT \
    -o $WDIR/mrf_seg_masked.nii.gz
SKIP

  ## Need to change label numbers to match segmentation protocol when computing dice
  if [[ $seg_mode == "coarse" ]]; then
  	$C3D_HOME/c3d $WDIR/mrf_seg_masked.nii.gz -replace 1 11 2 13 3 20 4 21 5 22 6 23 7 26 8 27 \
	9 37 10 42 11 45 12 46 13 51 14 53 15 59 16 5 17 61 18 6 -o $WDIR/mrf_seg_masked_relabel.nii.gz

  else
	  $C3D_HOME/c3d $WDIR/mrf_seg_masked.nii.gz -replace 1 11 2 12 3 13 4 17 5 20 6 21 7 22 8 23 \
        9 26 10 27 11 30 12 31 13 32 14 33 15 34 16 35 17 36 18 37 19 42 20 44 \
        21 45 22 46 23 4 24 51 25 53 26 56 27 57 28 59 29 5 30 61 \
        31 62 32 63 33 64 34 6 35 7 36 8 37 9  -o $WDIR/mrf_seg_masked_relabel.nii.gz
 
  fi

  # 4 5 6 7 8 9 11 12 13 17 20 21 22 23 26 27 30 31 32 33 34 35 36 37 42 44 45 46 51 53 56 57 59 61 62 63
 
  # Compute Dice with individual segmentations
  rm -rf $WDIR/gdsc.txt
  for seg in $(ls $WORK/histoseg_mapped_totemplate/*_histseg_to_template.nii.gz); do
    ID=$(basename $seg _histseg_to_template.nii.gz)

    $C3D_HOME/c3d $seg -as S -thresh 1 inf 1 0 \
      $WORK/vote_hfdb.nii.gz -thresh 1 1 1 0  -times -popas M \
      -push S $WDIR/mrf_seg_masked_relabel.nii.gz \
      -foreach -push M -times -endfor -label-overlap \
      > $WDIR/dice_report_$ID.txt

    echo "$WEIGHT $ID $(cat $WDIR/dice_report_$ID.txt| awk 'NR==3 {print $3}')" \
      >> $WDIR/gdsc.txt
  done

  # # Combine the subfields and the DB mask
  # $C3D_HOME/c3d $WDIR/mrf_seg_masked_relabel.nii.gz $WORK/vote_hfdb.nii.gz \
  #   -thresh 2 2 1 0 -as DB \
  #   -thresh 1 1 0 1 -times \
  #   -push DB -scale 5 -add \
  #   -o $WDIR/mrf_seg_withdb.nii.gz
}



function map_histo_segs_to_atlas()
{

  mode=${1?}

  seg_mode=${2?}

  N=$(cat $ROOT/scripts/subj_histo_annot.txt | wc -l)

  # Reslice individual segs
  mkdir -p $ROOT/histo_segs/dump
  for ((i=1;i<=$N;i++)); do
    fn=$(cat $ROOT/scripts/subj_histo_annot.txt | head -n $i | tail -n 1)
    subj=$(echo $fn | awk '{print $1}')
    $ROOT/scripts/job_launcher.sh -m 30G -o $ROOT/histo_segs/dump -N "hproc_${subj}" $0 process_seg $subj $mode $seg_mode
  done

<<"DONE"
  # Wait for completion
  $ROOT/scripts/job_launcher.sh -w "hproc_*" /bin/sleep 0

  # Combine the segs
  $ROOT/scripts/job_launcher.sh -m 40G -o $HDIR/dump -N "comb_${mode}" $0 combine_histo_subj $mode $seg_mode

  $ROOT/scripts/job_launcher.sh  -o $HDIR/dump -N "gdsc_${seg_mode}" $0 all_pair_histo_dice $seg_mode

  # Wait for completion
  $ROOT/scripts/job_launcher.sh -w "comb_${mode}" /bin/sleep 0

  if [[ $mode == "ncc" ]]; then

<<"FIN"
  # Determine best MRF weight
  for weight in $(echo 1 | awk '{ for (i=0.01; i<=0.101; i+=0.01) print i }'); do
	  $ROOT/scripts/job_launcher.sh -m 50G -o $HDIR/dump -N "mrf_${weight}" $ROOT/scripts/reslice_histo_segs.sh mrf ${weight} $seg_mode
  done

  $ROOT/scripts/job_launcher.sh -w "mrf_*" /bin/sleep 0
FIN

  # Link the chosen template
  CHOSEN_WEIGHT=010000
  
  ln -sf $TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}/mrf_vote/mrf_vote_${CHOSEN_WEIGHT}/mrf_seg_masked_relabel.nii.gz \
  $TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}/template_histo_seg_final.nii.gz

  fi

  ln -sf $TDIR/template_ncc/iter5/template_img.nii.gz $TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}/template_img.nii.gz

  # Create an image with only the hippocampus labels. Was used for visualizing the hippocampus surface only. 
  if [[ $seg_mode == "coarse" ]]; then
    c3d $TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}/template_histo_seg_final.nii.gz -as SEG \
    -thresh 1 inf 1 0 -as M -push SEG -retain-labels 5 13 14 15 20 23 21 22 26 37 42 44 46 \
    -thresh 1 inf 1 0 -push M \
    -add -o $TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}/template_hf_label.nii.gz

  fi

    # c3d -verbose $TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}/template_hf_label.nii.gz  \
    # -trim 5vox -replace 0 1000 \
    # -split \
    # -foreach -sdt -dup -times -sqrt -reciprocal -endfor \
    # -scale 0 -merge \
    # -replace 1000 0 -o $TEMPLATE_ROOTDIR/histology_annotation_${seg_mode}/seg_expanded.nii.gz
DONE


# The code below was used in Ravikumar et al. 2021 to plot the dispersion of the BA35/ERC boundary
<<"NOTNEEDED"
for ((i=1;i<=11;i++)); do
  	fn=$(cat $ROOT/scripts/subj_histo_annot.txt | head -n $i | tail -n 1)
 	subj=$(echo $fn | awk '{print $1}')

	# erc and ba35 boundary
 	c3d $HDIR/template_hist_29_coarse/totemplate/${subj}_histseg_to_template_phgmasked.nii.gz \
 	-as SEG -retain-labels 6 -replace 6 1 \
 	-dilate 1 1x1x1vox -as MASK \
 	-push SEG -retain-labels 27 \
 	-push MASK -times -thresh 27 27 1 0 -o $HDIR/template_hist_29_coarse/totemplate/${subj}_ba35_erc_boundary.nii.gz   
	# ba35 and ba36 boundary
        c3d $HDIR/template_hist_29_coarse/totemplate/${subj}_histseg_to_template_phgmasked.nii.gz \
        -as SEG -retain-labels 6 -replace 6 1 \
        -dilate 1 1x1x1vox -as MASK \
        -push SEG -retain-labels 11 \
        -push MASK -times -thresh 11 11 1 0 -o $HDIR/template_hist_29_coarse/totemplate/${subj}_ba35_ba36_boundary.nii.gz

         # ERC/paraSUB boundary
        c3d $HDIR/template_hist_29_coarse/totemplate/${subj}_histseg_to_template_phgmasked.nii.gz \
        -as SEG -retain-labels 27 -replace 27 1 \
        -dilate 1 1x1x1vox -as MASK \
        -push SEG -retain-labels 42 \
        -push MASK -times -thresh 42 42 1 0 -o $HDIR/template_hist_29_coarse/totemplate/${subj}_erc_parasub_boundary.nii.gz

done

SCALE=$(ls $HDIR/template_hist_29_coarse/totemplate/*_histseg_to_template_phgmasked.nii.gz | wc -l | awk '{print 1.0 / $1}')

#c3d $HDIR/template_hist_29_coarse/totemplate/*_ba35_erc_boundary.nii.gz \
#-accum -add -endaccum -scale $SCALE \
# -o $HDIR/template_hist_29_coarse/totemplate/boundary_dispersion_ba35_erc.nii.gz

#c3d $HDIR/template_hist_29_coarse/totemplate/*_ba35_ba36_boundary.nii.gz \
#-accum -add -endaccum -scale $SCALE \
# -o $HDIR/template_hist_29_coarse/totemplate/boundary_dispersion_ba35_ba36.nii.gz

c3d $HDIR/template_hist_29_coarse/totemplate/*_erc_parasub_boundary.nii.gz \
-accum -add -endaccum -scale $SCALE \
 -o $HDIR/template_hist_29_coarse/totemplate/boundary_dispersion_erc_parasub.nii.gz
NOTNEEDED

}


function main()
{

 # Run just once to copy relevant files from Paul's folder - updated code on 07/13/2022
 #fn_copy_paul_warps

 SUBJ_LIST=$ROOT/scripts/subj_histo_annot.txt
 #SUBJ_LIST=$ROOT/scripts/subj_map_hotspots.txt 
 
 N=$(cat $SUBJ_LIST | wc -l)
 mkdir -p $ROOT/histo_segs/dump

#coarse fine
<<"DONE"
for seg_mode in coarse fine; do
  for ((i=7;i<=7;i++)); do

    fn=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
    id=$(echo $fn | awk '{print $1}')

    #Loop through all blocks
    id_PY=$(cat $SUBJ_LIST | awk -v id=$id '{if ($1 == id) print $2}')
    read -r dummy blocks <<< "$(grep "^${id_PY}" "$ROOT/histo_segs/blockface_src.txt")"
      for block in $blocks; do
        echo $id $block
        #Only do this when downloading new scans from box
        #fn_copy_sydney_histo_annot $id $block

        # The functions below bridge the atlas and histology spaces
        if [[ -d  $ROOT/histo_segs/${id}/histo_annot/$block ]]; then
                
          #Map hotspot from atlas space to each subject's histology block for subjects where Sydney has completed
          # histology annotations
          #fn_map_hotspots_to_histo $id $block

          # Histo-To_MRI: Process entire histology segmentation per subject. Combine and interpolate
          $ROOT/scripts/job_launcher.sh -m 30G -o $HDIR/dump -N "map_to_raw_${id}_${block}" \
          $0 fn_map_histo_segs_to_hiresmri $id $block $seg_mode		

        fi            
              
              
      done
    done
done

    $ROOT/scripts/job_launcher.sh -w "map_to_raw_*" /bin/sleep 0

#Combine the segmentation across blocks 

for seg_mode in coarse fine; do
  for ((i=7;i<=7;i++)); do

      fn=$(cat $SUBJ_LIST | head -n $i | tail -n 1)
      id=$(echo $fn | awk '{print $1}')
      $ROOT/scripts/job_launcher.sh -m 30G -o $HDIR/dump -N "comb_${id}" \
          "$0" combine_blocks $id $seg_mode
  done 
done

# Wait for completion
$ROOT/scripts/job_launcher.sh -w "comb_*" /bin/sleep 0
DONE

for seg_mode in coarse fine; do
#Map the histology-based segmentations from raw to atlas space
  for mode in ncc; do # aff_ncc mst gshoot; do
    map_histo_segs_to_atlas $mode $seg_mode
  done
done


}


if [[ $1 ]]; then

  command=$1
  shift
        $command $@

else

  main
fi



