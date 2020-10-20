#!/bin/sh

inputpath="../pythia8244/dynamical_fragmentation/data/20200807_initial_condition_for_DCCI_pp_MB_CRon/multip_pp_pythia_weakStop.txt"
input_ext=""
outputdir="event_distrib_pp_pythia"
n=1

./analysis -n $n -output_dir ${outputdir} -input_path ${inputpath} #--ext ${input_ext}

