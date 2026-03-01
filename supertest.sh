#!/bin/bash

if [ $# != 3 ]; then
	echo "Usage: supertest.sh <a.out> <r> <n(до какого n делать тесты)>"
else
    exe=$1
    r=$2
    n=$3
    
    if [ -f $exe ];then
    echo
    else
        echo "$exe doesn't exist"
        exit
    fi

    for ((i=2; i<=$n; i++)); do

        for ((j=1; j<i; j++)); do

            for ((k=1; k<4;k++)); do

		for ((q = 1; q <=4 ; q++)); do
		    echo "= mpirun -np $q $exe $i $j $r $k ="	
                    echo "=============================== N = $i ==== M = $j ==== P = $q ==== S = $k ========================"
                    mpirun -np $q $exe $i $j $r $k
                    echo "---------------------------------------------------------------------------------------------------"
		    echo -e "\n\n"
		done

            done

        done

    done

fi
