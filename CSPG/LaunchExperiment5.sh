#!/bin/bash
#title           :LaunchExperiment1.sh
#description     :This script tests the influence of the first approximation level.
#author		 :Jean-Luc Bouchot
#date            :2019/06/05
#version         :0.1    
#usage		 :bash LaunchExperiment1.sh
#notes           :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Script to run the fifth experiment
# This experiment tests the effect of the sparsity level picked on the solution
#
# We keep all the values fixed except for the constant appearing in s_L.

 
# List of parameters
d=20 # number of cosines
start_h0=10 # will be used as the discretization step for the first level. 
Lmax=6 # Target number of discretization steps -> Let this vary! 
J=3
algo=wiht # Which algorithm should be used
vj=1.08 # Value of the constant coefficients
nbSamples=new # What should be the number of samples 
nbtests=100 # A few tests at the end to make sure it somewhat worked
powerTrig=4.5 # Power of the trigonometric decay
abar=10 # Constant mean field
flucImportance=1 # Importance of the fluctuations
dotensor=TRUE # Use a tensor-based computation instead of building the whole sensing matrix
wCosine=0.25 # How important the 'j' component in the decay is
p0=0.25 # Compressibility in the original space
p=0.3 # Compressibility in the smoothness scale
sJ=40 # Constant used for the first level of approximation



expBasename="Exp5Dim2WCosine${d}InfluenceSparsity"


sLmax=80 # Constant appearing in front of the number of samples of the target discretization
sL=15

while (( $sL <= $sLmax ))
do
	echo "RUNNING THE EXPERIMENT WITH J = $J"
	folder=$expBasename$J
	python test_wCosine_2D_avg_v_ML.py -d $d -o WeightedCosine2D -L $Lmax -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder
	sL=$(( $sL + 3 ))
done
