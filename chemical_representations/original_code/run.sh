#!/bin/bash

max=2
count=0

while [[ $count -lt 22 ]]
do
    pnum=0
    for i in $(pgrep -U liz python); do
        pnum=$(( $pnum + 1 ))
    done

    if [[ $pnum -lt $max ]]
    then
        python DM_v_1_generation.py -f DM_data/source_files/chemical_ids_$count.txt&
        echo /DM_data/source_files/chemical_ids_$count.txt
        count=$(($count+1))
    fi;
    sleep 5m
done
