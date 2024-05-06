#!/bin/bash
#$ -S /bin/bash

ROOT_DIR=/project/hippogang_3/sravikumar
ADIR_ROOT=$ROOT_DIR/atlasPHG2019

# PATH to CMREP commands
#cmrep version with meshglm with nuissance variables (fl) and missing data
# and binary mesh_merge_array
CMREP_HOME=$ROOT_DIR/packages/cmrepgccdbg_withfl/

# Cmrep build with lmshoot with PLY compatability
CMREP_LMSHOOT=$ROOT_DIR/atlasPHG2019/pkgs/cmrep-1.0.0-Linux-x86_64_nov2021/bin
export PATH=$CMREP_LMSHOOT:$PATH

export LD_LIBRARY_PATH=$CMREP_HOME/lib:$HOME/lib:$LD_LIBRARY_PATH

# PATH to C3D
C3D_HOME=$ADIR_ROOT/pkgs

# PATH to Greedy
#GREEDY_HOME=/data/picsl-build/pauly/greedy/gcc64rel
#GREEDY_HOME=$ADIR/pkgs
GREEDY_HOME=$ADIR_ROOT/pkgs/greedy-1.2.0-Linux-gcc64/bin

module load ANTs
module load R

#Activate python environment with pymeshlab etc. 
module load python/3.6.1
source $ADIR_ROOT/pkgs/venv_python/bin/activate

# Path to ANTS
#ANTSDIR=/share/apps/ANTs/2014-06-23/build/bin

# Path to R
#RHOME=/share/apps/R/R-3.2.3/

MESH_SRCODE=$ADIR_ROOT/scripts/mesh_fill_missing_data/bin

MESH2IMG=$ADIR_ROOT/scripts/mesh_to_image/bin

# Make sure there is a tmpdir
if [[ ! $TMPDIR ]]; then
  TMPDIR=/tmp/atlas_${PPID}
  mkdir -p $TMPDIR
fi

