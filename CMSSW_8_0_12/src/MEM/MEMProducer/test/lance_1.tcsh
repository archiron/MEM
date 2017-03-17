#!/bin/bash

DATE=$(date +'%Y.%m.%d_%H.%M.%S')
echo $DATE

mkdir $DATE
cd $DATE

nb_evts=200
if [ "$1" == "" ] 
then
	nb_evts=0
else
    nb_evts=$1
    echo "nb evts = $nb_evts"
fi

#mpirun -np 1 cmsRun ../MEMProducer.write_cfg.py $nb_evts 
cmsRun ../MEMProducer.write_cfg.py $nb_evts 

echo "$DATE fini !"
