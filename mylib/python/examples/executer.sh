#!/bin/bash

if [ "$1" = "help" ]; then
echo "This script must be run as follows:"
echo "sh executer.sh reactor problem mode"
exit 0
fi

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
    gdb --args ../../../testing/main.exe -f ${reactor}/${problem}
elif [ "${mode}" = "memory" ]; then
    echo "Running and valgrind"
    valgrind --leak-check=yes ../../../testing/main.exe -f ${reactor}/${problem}
else 
    echo "Running"
    ../../../testing/main.exe -f ${reactor}/${problem}
fi

