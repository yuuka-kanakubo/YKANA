#!/bin/sh

inputpath="../pythia8244/dynamical_fragmentation/data/20200903_initial_condition_for_DCCI_PbPb_MB_CRon_pTRef1.0/multip_PbPb_pT0Ref1.0.txt"
input_ext=""
outputdir="event_distrib_PbPb_pT0Ref1.0_pythia"
n=1

./analysis -n $n -output_dir ${outputdir} -input_path ${inputpath} #--ext ${input_ext}

