#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# NB At the moment you have to manually make sure you only submit the correct number of jobs for number of processors
#
# This script assumes that the following has been run successfully:
# scons compile_only=1 build=GccOpt test_suite=projects/SaraJaneD/test/TestSimulationRunsVaryingTDSpringStrength.hpp
#

TRANSIT_DIFFERENTIATED_MULTIPLIER[0]="1.0"
TRANSIT_DIFFERENTIATED_MULTIPLIER[1]="1.1"
TRANSIT_DIFFERENTIATED_MULTIPLIER[2]="1.2"
TRANSIT_DIFFERENTIATED_MULTIPLIER[3]="1.3"
TRANSIT_DIFFERENTIATED_MULTIPLIER[4]="1.4"
TRANSIT_DIFFERENTIATED_MULTIPLIER[5]="1.5"
TRANSIT_DIFFERENTIATED_MULTIPLIER[6]="1.6"
TRANSIT_DIFFERENTIATED_MULTIPLIER[7]="1.7"
TRANSIT_DIFFERENTIATED_MULTIPLIER[8]="1.8"
TRANSIT_DIFFERENTIATED_MULTIPLIER[9]="1.9"
TRANSIT_DIFFERENTIATED_MULTIPLIER[10]="2.0"
TRANSIT_DIFFERENTIATED_MULTIPLIER[11]="2.1"
TRANSIT_DIFFERENTIATED_MULTIPLIER[12]="2.2"
TRANSIT_DIFFERENTIATED_MULTIPLIER[13]="2.3"
TRANSIT_DIFFERENTIATED_MULTIPLIER[14]="2.4"
TRANSIT_DIFFERENTIATED_MULTIPLIER[15]="2.5"
TRANSIT_DIFFERENTIATED_MULTIPLIER[16]="2.6"
TRANSIT_DIFFERENTIATED_MULTIPLIER[17]="2.7"
TRANSIT_DIFFERENTIATED_MULTIPLIER[18]="2.8"
TRANSIT_DIFFERENTIATED_MULTIPLIER[19]="2.9"
TRANSIT_DIFFERENTIATED_MULTIPLIER[20]="3.0"
TRANSIT_DIFFERENTIATED_MULTIPLIER[21]="3.1"
TRANSIT_DIFFERENTIATED_MULTIPLIER[22]="3.2"
TRANSIT_DIFFERENTIATED_MULTIPLIER[23]="3.3"
TRANSIT_DIFFERENTIATED_MULTIPLIER[24]="3.4"
TRANSIT_DIFFERENTIATED_MULTIPLIER[25]="3.5"
TRANSIT_DIFFERENTIATED_MULTIPLIER[26]="3.6"
TRANSIT_DIFFERENTIATED_MULTIPLIER[27]="3.7"
TRANSIT_DIFFERENTIATED_MULTIPLIER[28]="3.8"
TRANSIT_DIFFERENTIATED_MULTIPLIER[29]="3.9"
TRANSIT_DIFFERENTIATED_MULTIPLIER[30]="4.0"

for (( i=0 ; i<${#TRANSIT_DIFFERENTIATED_MULTIPLIER[*]} ; i++))
do
	echo "Beginning run for transit-differentiated spring strength = ${TRANSIT_DIFFERENTIATED_MULTIPLIER[$i]}."
    # NB "nice -20" gives the jobs low priority (good if they are going to dominate and and no slower if nothing else is going on)
    # ">" directs std::cout to the file.
    # "2>&1" directs std::cerr to the same place.
    # "&" on the end lets the script carry on and not wait until this has finished.
   	nice -20 ../build/optimised/TestSimulationRunsVaryingTDSpringStrengthRunner ${TRANSIT_DIFFERENTIATED_MULTIPLIER[$i]} > box_model_${TRANSIT_DIFFERENTIATED_MULTIPLIER[$i]}_output.txt 2>&1 &
done

echo "Jobs submitted"