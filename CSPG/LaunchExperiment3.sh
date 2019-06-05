#!/bin/bash

# Script to run the third experiment
# This experiment tests the importance of the target discretization level
#
# For a two dimensional (spatial) problem, we will set the number of discretization levels
# and the original grid size and let the target approximation level run from 3 to Lmax, 
# while the first level discretization will follow J = L-3

 
# List of parameters
d=20 # number of cosines
start_h0=10 # will be used as the discretization step for the first level. 
Lmax=4 # Target number of discretization steps 
Lmin=3 # minimum number of multi-levels computed
algo=wiht # Which algorithm should be used
vj=1.08 # Value of the constant coefficients
nbSamples=new # What should be the number of samples 
nbtests=100 # A few tests at the end to make sure it somewhat worked
powerTrig=4.5 # Power of the trigonometric decay
abar=10 # Constant mean field
flucImportance=1 # Importance of the fluctuations
sL=40 # Constant appearing in front of the number of samples of the target discretization
dotensor=TRUE # Use a tensor-based computation instead of building the whole sensing matrix
wCosine=0.25 # How important the 'j' component in the decay is
p0=0.25 # Compressibility in the original space
p=0.3 # Compressibility in the smoothness scale
sJ=40 # Constant used for the first level of approximation

expBasename=Exp3Dim2WCosine${d}InfluenceTarget

for ((L=$Lmin; L<=$Lmax; L++))
do
	echo "RUNNING THE EXPERIMENT WITH L = $L"
	folder=$expBasename$L
        J=$(($L-$Lmin))
        echo "Starting index is J = $J with end index L = $L" 
	python test_wCosine_2D_avg_v_ML.py -d $d -o WeightedCosine2D -L $L -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder
done
