#!/bin/bash
# Try to run the basic case of the model
# give a run name as an input

export CASE=$1
export CASEROOT=$HOME/projects/def-czg/bps5238/cesmruns/$CASE
#export CASEROOT=$HOME/scratch/cesmruns/$CASE
export MACH=computecanada
export RES=f19_g16
#export RES=T62_g16
export COMPSET=E_1850_CAM5
#export COMPSET=X
#export COMPSET=E_1850_CAM5
#export COMPSET=D_NORMAL_YEAR
module load cesm\/1_2_2 
create_newcase --case $CASEROOT --res $RES --compset $COMPSET --mach $MACH








