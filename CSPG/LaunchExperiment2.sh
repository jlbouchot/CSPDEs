#!/bin/bash
#title           :LaunchExperiment2.sh
#description     :This script tests the influence of the target approximation level.
#author		 :Jean-Luc Bouchot
#date            :2019/06/05
#version         :0.1    
#usage		 :bash LaunchExperiment2.sh
#notes           :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Script to run the second experiment
# This experiment tests the importance of the target discretization level
#
# For a two dimensional (spatial) problem, we will set the first discretization level
# and the original grid size and let the target approximation level run from J to Lmax.

 
# List of parameters
d=20 # number of cosines
start_h0=40 # will be used as the discretization step for the first level. 
Lmax=8 # Target number of discretization steps -> Let this vary! 
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
J=0 # Where to do the first (single level) approximation

expBasename=Exp2Dim2WCosine${d}InfluenceLmax

for ((L=$J; L<=$Lmax; L++))
do
	echo "RUNNING THE EXPERIMENT WITH L = $L"
	folder=$expBasename$L
	python test_wCosine_2D_avg_v_ML.py -d $d -o WeightedCosine2D -L $L -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder
done
