#!/bin/bash

if [ $# != 2 ]; then
    echo "Usage: <exec> <out file>"

else
    exe=$1
    outf=$2

    if [ -f $outf ]; then
        echo -n
    else
        echo "File $outf doesn't exist"
        echo "File $outf created"
        touch $outf
    fi

    if [ -f $exe ]; then
        echo -n
    else
        echo "Executable $exe doesn't exist"
        exit 0
    fi

    greptest=$(grep "TEST" $outf)
    # echo $greptest

    if [[ $greptest =~ "TEST" ]]; then
        numtest="${greptest: -1:1}"
        numtest=$((numtest+1))
        # echo $numtest
        echo "************************************************************************************************" >> $outf
        echo "************************************************************************************************" >> $outf
        echo
        echo TEST $numtest >> $outf
    else
        numtest="0";
        echo TEST $numtest > $outf
    fi

    for((i=9000; i<=15000; i+=2000)); do
        echo "= mpirun -np 4 $exe $i 60 5 1 ="	>> $outf
        echo "=============================== N = $i ==== M = 60 ==== P = 4 ==== S = 1 ========================" >> $outf
        mpirun -np 4 $exe $i 60 5 1 >> $outf
        echo "----------------------------------------------------------------------------------" >> $outf
    done

     
    

fi