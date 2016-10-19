#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#

start_run=1; #starting index
end_run=5; #final index

for ((i=start_run; i<=end_run; i++))
do
echo "Beginning runs from ${start_run} to ${end_run}."
# NB "nice -20" gives the jobs low priority (good if they are going to dominate and and no slower if nothing else is going on)
# ">" directs std::cout to the file.
# "2>&1" directs std::cerr to the same place.
# "&" on the end lets the script carry on and not wait until this has finished.
nice ../build/optimised/TestOrganoid/TestOrganoidRingWithPanethCellsRunner -myintval ${i} &
done

echo "Jobs submitted"