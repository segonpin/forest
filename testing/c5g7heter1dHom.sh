#!/bin/bash

reactor="c5g7heter1d"
mode=$1

problem="inputHomAss"
sh executer.sh ${reactor} ${problem} ${mode}

problem="inputHomPin"
sh executer.sh ${reactor} ${problem} ${mode}
