#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# NB At the moment you have to manually make sure you only submit the correct number of jobs for number of processors
#
# This script assumes that the following has been run successfully:
# scons compile_only=1 build=GccOpt test_suite=projects/SaraJaneD/test/TestBoxModel/TestSimulationRunsForBoxModel.hpp
#
BASAL_LAMINA_PARAMETER[0]="12.0"
#BASAL_LAMINA_PARAMETER[1]="8.0"
#BASAL_LAMINA_PARAMETER[2]="2.0"
#BASAL_LAMINA_PARAMETER[3]="3.0"
#BASAL_LAMINA_PARAMETER[4]="4.0"
#BASAL_LAMINA_PARAMETER[5]="5.0"
#BASAL_LAMINA_PARAMETER[6]="6.0"
#BASAL_LAMINA_PARAMETER[7]="7.0"
#BASAL_LAMINA_PARAMETER[8]="8.0"
#BASAL_LAMINA_PARAMETER[9]="9.0"
#BASAL_LAMINA_PARAMETER[10]="10.0"
#BASAL_LAMINA_PARAMETER[11]="11.0"
#BASAL_LAMINA_PARAMETER[12]="12.0"
#BASAL_LAMINA_PARAMETER[13]="13.0"
#BASAL_LAMINA_PARAMETER[14]="14.0"
#BASAL_LAMINA_PARAMETER[15]="15.0"
#BASAL_LAMINA_PARAMETER[16]="16.0"
#BASAL_LAMINA_PARAMETER[17]="17.0"
#BASAL_LAMINA_PARAMETER[18]="18.0"
#BASAL_LAMINA_PARAMETER[19]="19.0"
#BASAL_LAMINA_PARAMETER[20]="20.0"

for (( i=0 ; i<${BASAL_LAMINA_PARAMETER[*]} ; i++))
do
	echo "Beginning run for basal lamina parameter = ${BASAL_LAMINA_PARAMETER[$i]}."
    # NB "nice -20" gives the jobs low priority (good if they are going to dominate and and no slower if nothing else is going on)
    # ">" directs std::cout to the file.
    # "2>&1" directs std::cerr to the same place.
    # "&" on the end lets the script carry on and not wait until this has finished.
   	nice -20 ../build/optimised/TestBoxModel/TestSimulationRunsForBoxModelRunner ${BASAL_LAMINA_PARAMETER[$i]} > box_model_${BASAL_LAMINA_PARAMETER[$i]}_output.txt 2>&1 &
done

echo "Jobs submitted"
