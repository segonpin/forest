#!/bin/bash

echo 'GeorgiaTech1D assemblies'

reactor="GeorgiaTech1D"
mode=$1

problem="GeorgiaTech1DAss1"
sh executer.sh ${reactor} ${problem} ${mode}

problem="GeorgiaTech1DAss2"
sh executer.sh ${reactor} ${problem} ${mode}

problem="GeorgiaTech1DAss3"
sh executer.sh ${reactor} ${problem} ${mode}

problem="GeorgiaTech1DAss4"
sh executer.sh ${reactor} ${problem} ${mode}



