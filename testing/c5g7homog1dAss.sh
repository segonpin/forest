#!/bin/bash

echo 'C5G7 homog 1D assemblies'

reactor="c5g7homog1d"
mode=$1

problem="inputAss1"
sh executer.sh ${reactor} ${problem} ${mode}

problem="inputAss2"
sh executer.sh ${reactor} ${problem} ${mode}



