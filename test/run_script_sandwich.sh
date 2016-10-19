#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# NB At the moment you have to manually make sure you only submit the correct number of jobs for number of processors
#
# This script assumes that the following has been run successfully:
# scons compile_only=1 build=GccOpt test_suite=projects/SaraJaneD/test/Test2DCrossSectionModel/TestMultipleSimulationCrossSectionModel.hpp
#

BASAL_LAMINA_PARAMETER[0]="12.0"
#BASAL_LAMINA_PARAMETER[1]="13.0"
#BASAL_LAMINA_PARAMETER[2]="14.0"
#BASAL_LAMINA_PARAMETER[3]="15.0"
#BASAL_LAMINA_PARAMETER[4]="16.0"
#BASAL_LAMINA_PARAMETER[5]="17.0"
#BASAL_LAMINA_PARAMETER[6]="18.0"
#BASAL_LAMINA_PARAMETER[7]="19.0"
#BASAL_LAMINA_PARAMETER[8]="20.0"

for (( i=0 ; i<${#BASAL_LAMINA_PARAMETER[*]} ; i++))
do
	echo "Beginning run for basal lamina parameter = ${BASAL_LAMINA_PARAMETER[$i]}."
    # NB "nice -20" gives the jobs low priority (good if they are going to dominate and and no slower if nothing else is going on)
    # ">" directs std::cout to the file.
    # "2>&1" directs std::cerr to the same place.
    # "&" on the end lets the script carry on and not wait until this has finished.
   	nice -20 ../build/optimised/Test2DCrossSectionModel/TestMultipleSimulationsCrossSectionModelRunner ${BASAL_LAMINA_PARAMETER[$i]} > box_model_${BASAL_LAMINA_PARAMETER[$i]}_output.txt 2>&1 &
done

echo "Jobs submitted"
