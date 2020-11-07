#!/bin/sh

inputpath="../pythia8244/dynamical_fragmentation/data/20200917_initial_condition_for_DCCI_pp_MB_CRon_pTRef1.4/multip_pp_pythia_pTRef1.4_weakStop.txt"
input_ext=""
outputdir="event_distrib_pp_pythia_pT0Ref1.4"
n=1

./analysis -n $n -output_dir ${outputdir} -input_path ${inputpath} #--ext ${input_ext}

