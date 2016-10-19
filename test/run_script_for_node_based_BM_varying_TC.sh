#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#

#start_run=6; #starting index
#end_run=10; #final index

#Define the values of beta to run in the simulation
target_curvature_values[0]="0.0"
target_curvature_values[1]="0.1"
target_curvature_values[2]="0.2"
target_curvature_values[3]="0.3"
target_curvature_values[4]="0.4"


for ((i=0; i<${#beta_values[*]}; i++))
do
echo "Running test for 1/R = ${target_curvature_values[$i]}."
# NB "nice -20" gives the jobs low priority (good if they are going to dominate and and no slower if nothing else is going on)
# ">" directs std::cout to the file.
# "2>&1" directs std::cerr to the same place.
# "&" on the end lets the script carry on and not wait until this has finished.
nice ../build/optimised/TestOrganoid/TestNodeBasedBasementMembraneModelRunner -mydoubleval $target_curvature_values[$i]} > node_based_box_model_TC_${taget_curvature_values[$i]}_output.txt 2>&1 &
done

echo "Job submitted"
