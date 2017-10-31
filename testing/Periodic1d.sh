#!/bin/bash

reactor="Periodic1D"

# Heterogeneous

problem="PeriodicHeterogeneous1D"
mode=$1
sh executer.sh ${reactor} ${problem} ${mode}

problem="PeriodicHeterogeneous1DAss"
mode=$1
sh executer.sh ${reactor} ${problem} ${mode}

# Homogeneous

problem="PeriodicHomogeneous1D"
mode=$1
sh executer.sh ${reactor} ${problem} ${mode}

problem="PeriodicHomogeneous1DAss"
mode=$1
sh executer.sh ${reactor} ${problem} ${mode}

