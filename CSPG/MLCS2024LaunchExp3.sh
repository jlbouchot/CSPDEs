#!/bin/bash
#title           :LaunchExperiment4.sh
#description     :This script tests the influence of the dimensionality of the parameters.
#author		 :Jean-Luc Bouchot
#date            :2019/07/09
#version         :0.1    
#usage		 :bash LaunchExperiment4.sh
#notes           :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# Script to run the third experiment
# This experiment tests the importance of the target discretization level
#
# For a two dimensional (spatial) problem, we will set the number of discretization levels
# and the original grid size and let the target approximation level run from 3 to Lmax, 
# while the first level discretization will follow J = L-3

 
# List of parameters
# dmax=30 # number of cosines
start_h0=10 # will be used as the discretization step for the first level. 
Lmax=6 # Target number of discretization steps 
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

expBasename=Exp4Dim2WCosineDimensionalityd

for d in 8 10 13 16 20 25 32
#for ((d=$Lmin; L<=$Lmax; L++))
do
	echo "RUNNING THE EXPERIMENT WITH d = $d"
	folder=$expBasename$d
	python test_wCosine_2D_avg_v_ML.py -d $d -o WeightedCosine2D -L $Lmax -s $Lmin -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder
done
