#!/bin/bash

mode=$1

#echo "c5g7homog1d"
#sh c5g7homog1d.sh ${mode}
#sh c5g7homog1dAss.sh ${mode}
#sh c5g7homog1dHom.sh ${mode}

echo "c5g7heter1d"
sh c5g7heter1d.sh ${mode}
sh c5g7heter1dAss.sh ${mode}
sh c5g7heter1dHom.sh ${mode}

