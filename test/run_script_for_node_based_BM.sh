#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#

#start_run=6; #starting index
#end_run=10; #final index

#Define the values of beta to run in the simulation
beta_values[0]="0"
beta_values[1]="1.0"
beta_values[2]="2.0"
beta_values[3]="3.0"
beta_values[4]="4.0"
beta_values[5]="5.0"
beta_values[6]="6.0"
beta_values[7]="7.0"
beta_values[8]="8.0"
beta_values[9]="9.0"
beta_values[10]="10.0"

for ((i=0; i<${#beta_values[*]}; i++)) 
do
echo "Running test for Beta = ${beta_values[$i]}."
# NB "nice -20" gives the jobs low priority (good if they are going to dominate and and no slower if nothing else is going on)
# ">" directs std::cout to the file.
# "2>&1" directs std::cerr to the same place.
# "&" on the end lets the script carry on and not wait until this has finished.
nice ../build/optimised/TestOrganoid/TestNodeBasedBasementMembraneModelRunner -mydoubleval ${beta_values[$i]} > node_based_box_model_${beta_values[$i]}_output.txt 2>&1 &
done

echo "Job submitted"
