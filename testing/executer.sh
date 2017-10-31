#!/bin/bash

reactor=$1
problem=$2
if [ "$3" = "debug" ]; then
mode="debug"
elif [ "$3" = "memory" ]; then
mode="memory"
else 
mode="running"
fi

echo "REACTOR: ${reactor}; PROBLEM: ${problem}; MODE: ${mode}"

if [ "${mode}" = "debug" ]; then
    echo "Running and debuging"
    gdb --args ./main.exe -f ../reactors/${reactor}/${problem}
elif [ "${mode}" = "memory" ]; then
    echo "Running and valgrind"
    valgrind --leak-check=yes ./main.exe -f ../reactors/${reactor}/${problem}
else 
    echo "Running"
    ./main.exe -f ../reactors/${reactor}/${problem}
fi

