#!/bin/bash

# Script to run the fourth experiment
# This experiment tests the influence of the fluctuations on the parametric solution
#
# For a two dimensional (spatial) problem, we will set all parameters but the relative importance of the fluctuations

 
# List of parameters
d=20 # number of cosines
start_h0=10 # will be used as the discretization step for the first level. 
Lmax=7 # Target number of discretization steps 
Lmin=1 # minimum number of multi-levels computed
algo=wiht # Which algorithm should be used
vj=1.08 # Value of the constant coefficients
nbSamples=new # What should be the number of samples 
nbtests=100 # A few tests at the end to make sure it somewhat worked
powerTrig=4.5 # Power of the trigonometric decay
abar=10 # Constant mean field
sL=40 # Constant appearing in front of the number of samples of the target discretization
dotensor=TRUE # Use a tensor-based computation instead of building the whole sensing matrix
wCosine=0.25 # How important the 'j' component in the decay is
p0=0.25 # Compressibility in the original space
p=0.3 # Compressibility in the smoothness scale
sJ=40 # Constant used for the first level of approximation

# How are we going to deal with the various importances? 
base=10
imin=1
imax=15

expBasename=Exp3Dim2WCosine${d}InfluenceFluctuations

for ((curImp=$imin; curImp<=$imax; curImp++))
do
	impValue=$(echo "scale=2; $curImp/$base" | bc)
	echo "RUNNING THE EXPERIMENT WITH IMPORTANCE = $impValue"
	folder=$expBasename$curImp
        echo "Starting index is J = $J with end index L = $L" 
	python test_wCosine_2D_avg_v_ML.py -d $d -o WeightedCosine2D -L $Lmax -s $Lmin -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $impValue -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder
done
