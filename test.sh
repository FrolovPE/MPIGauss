#!/bin/bash

# echo "Usage: test.sh <executable(./a.out, './' before file must have)> <dir(example matr/ , '/' after dir name must have )>"
# sleep 3


if [ $# != 2 ]; then
	echo "Usage: test.sh <executable(./a.out or a.out)> <dir(example matr/ , '/' is must have )>"
else
	
	 exe=$1
	 dir=$2

	

	 if [ -f $exe ];then
		echo -n
	else 
		echo "executable file doesn't exist"
		exit
	fi

	if [ "./" == *"$exe"* ]; then
	 	echo -n
	else
		exe="./$exe"
	fi

	if [ -d $dir ];then
		echo -n
	else 
		echo "dir doesn't exist"
		exit
	fi


	if [ "/" == *"$dir"* ]; then
	 	echo -n
	else
		dir="$dir/"
	fi

	 for i in $dir*; do
	     ii="${i##*/}"
	     filename="${i%.*}"
	     filename="${filename##*/}"
	     ext="${i##*.}"
#	     echo "$ext"
#	     echo $i
		if [ $ext != "sh" ] && [ $ii != "3" ] && [ $ext != "." ] && [ $ext != "zip" ]; then
		         if [ $filename == "c" ] || [ $filename == "d" ] || [ $filename == "e" ]; then
				 echo $i
				 		
						for((k = 1; k <= 6 ; k++)); do
							for((q = 1 ; q <=4; q++)); do
						echo "------------------------------------------------------------------------"
						echo "mpirun -np $q $exe 6 $k 6 0 $i"
		                 mpirun -np $q $exe 6 $k 6 0 $i
		                 echo ""
							done
						 done
		         else
		             echo ""
		             echo $i
						for((k = 1; k <= 4 ; k++)); do
                                                        for((q = 1 ; q <=4; q++)); do
                                                echo "------------------------------------------------------------------------"
                                		 	echo "mpirun -np $q $exe 4 $k 4 0 $i"
                                			 mpirun -np $q $exe 4 $k 4 0 $i
                                		 	echo ""
                		                        done
                                                 done
		         fi
		fi
		echo
	 done
fi
