#!/bin/bash
#title           :MLCSPG2024LaunchExp4.sh
#description     :This script tests the influence of the target approximation level.
#author		 :Jean-Luc Bouchot
#date            :2024/09/05
#version         :0.1    
#usage		 :bash MLCSPG2024LaunchExp4.sh
#notes           :Install FEniCS, CVXPY, progressbar before using.
#==============================================================================

# This experiment tests the importance of the target discretization level
#
# For a two dimensional (spatial) problem, we will set the number of discretization levels
# and the original grid size and let the target approximation level run from 2 to Lmax, 
# while the first level discretization will follow J = L-2

 
# List of parameters
d=10 # number of cosines
start_h0=20 # will be used as the discretization step for the first level. 
L=6 # Target number of discretization steps 
Jmin=0 # minimum number of multi-levels computed
Jmax=3
algo=wiht # Which algorithm should be used
vj=1.005 # Value of the constant coefficients
nbSamples=new # What should be the number of samples 
nbtests=50 # A few tests at the end to make sure it somewhat worked
powerTrig=5 #4 # Power of the trigonometric decay
abar=10 # Constant mean field
flucImportance=1 # Importance of the fluctuations
sL=80 # Constant appearing in front of the number of samples of the target discretization
dotensor=True # Use a tensor-based computation instead of building the whole sensing matrix
wCosine=1 # How important the 'j' component in the decay is
p0=0.33 #0.5 # Compressibility in the original space
p=0.33 #0.5 # Compressibility in the smoothness scale
sJ=5 # Constant used for the first level of approximation
exponent=0.166 #0.25
expBasename=Exp4H0${start_h0}Dim2WCosine${d}InfluenceTarget${L}J
nbDetails=2

for ((J=$Jmin; J<=$Jmax; J++))
do
	echo "RUNNING THE EXPERIMENT WITH J = $J"
	folder=$expBasename$J
	mkdir -p $folder
	#J=$(($L-$nbDetails))
	echo "Starting index is J = $J with end index L = $L" 
	# python test_wCosine_2D_avg_p_ML.py -d $d -o WeightedCosine2D -L $L -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder -E $exponent -k False > $folder/stdoutput.txt 
	python test_wCosine_2D_avg_p_ML.py -d $d -o WeightedCosine2D -L $L -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder -E $exponent > $folder/stdoutput.txt # -k False
	# python test_wCosine_2D_avg_p_ML.py -d $d -o WeightedCosine2D -L $L -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder -E $exponent > $folder/stdoutput.txt # -k False
	# python test_wCosine_2D_avg_p_ML.py -d $d -o WeightedCosine2D -L $L -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder -E $exponent -k False 
	# python test_wCosine_2D_avg_p_ML.py -d $d -o WeightedCosine2D -L $L -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder
done
