#!/bin/sh

inputpath="../pythia8244/default_pythia_mymain/data/20201020_default_pPb_WeakDecStop/centrality_cut/multip_pPb_pythia_weakStop.txt"
input_ext=""
outputdir="event_distrib_pp_default_mymain"
n=1

./analysis -n $n -output_dir ${outputdir} -input_path ${inputpath} #--ext ${input_ext}

