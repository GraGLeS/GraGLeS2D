if [ $# -lt 4 ]; then
	echo "Usage: BenchmarkJobs.sh <executable_filename> <benchmark_parametes_filename> <Email> <NumThreads1,MaxMem,TimeLimit> [NumThreads2] ... [NumThreadsN]"
	exit 2
fi

if ! [ -x $1 ]; then
	echo "$1 is not an executable file!"
	exit 2
fi

if ! [ -f $2 ]; then
	echo "$2 was not found or is not a regular file!"
	exit 2
fi

if [ -z $3 ]; then
	echo "Specify an e-mail!"
	exit 2
fi

EXECUTABLE_FILENAME=$1
PARAMETERS_FILENAME=$2
EMAIL=$3
JOBSCRIPT_NAME_TEMPLATE="JobScript_Threads"

for i in "${@:4}"; do
	IFS=',' read -ra arr <<< "$i"    
    if [ -z ${arr[0]} ]; then
    	echo $i" has an invalid thread count"
    fi
    
    if [ -z ${arr[1]} ]; then
    	echo $i" has an invalid max memory"
    fi
    
    if [ -z ${arr[2]} ]; then
    	echo $i" has an invalid time limit"
    fi
       
    JOBSCRIPT_NAME=$JOBSCRIPT_NAME_TEMPLATE'_'${arr[0]}.sh
    mkdir Simulation_${arr[0]}_Threads
    cp $EXECUTABLE_FILENAME Simulation_${arr[0]}_Threads
    cp $PARAMETERS_FILENAME Simulation_${arr[0]}_Threads
    cp lib* Simulation_${arr[0]}_Threads
    cd Simulation_${arr[0]}_Threads
 
    echo "#!/usr/bin/env zsh" > $JOBSCRIPT_NAME
    echo "#BSUB -J LEVELSET" >> $JOBSCRIPT_NAME
    echo "#BSUB -o LEVELSET.%J" >> $JOBSCRIPT_NAME
    echo "#BSUB -e LEVELSET.e%J" >> $JOBSCRIPT_NAME
    echo "#BSUB -M ${arr[1]}" >> $JOBSCRIPT_NAME
    echo "#BSUB -W ${arr[2]}" >> $JOBSCRIPT_NAME
    echo "#BSUB -u "$EMAIL >> $JOBSCRIPT_NAME
    echo "#BSUB -N" >> $JOBSCRIPT_NAME
    echo "#BSUB -n "${arr[0]} >> $JOBSCRIPT_NAME
    echo "#BSUB -a \"bcs openmp\" ">> $JOBSCRIPT_NAME
    echo "./$EXECUTABLE_FILENAME $PARAMETERS_FILENAME" >> $JOBSCRIPT_NAME
    
    bsub < $JOBSCRIPT_NAME
    cd ..
    #rm $JOBSCRIPT_NAME
done
