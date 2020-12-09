#!/bin/sh


inputpath="../DCCI/data/20201204_PP_PT0REF2.0_SIGMA0.5"
input_ext="/hadronFinal_corecorona_weakStop.txt"
outputdir="debug"
n=10000

#Do not modify this.
#---------------------
log_dname="log/"
data_dir="data/"
today=$(date "+%Y%m%d")
log_fname=${today}${dirname}
if [ ! -d $log_dname ]
then
    echo "Directory "$log_dname" DOES NOT exists." 
    echo "mkdir "$log_dname
    mkdir $log_dname
fi
if [ ! -d $data_dir ]
then
    echo "Directory "$data_dir" DOES NOT exists." 
    mkdir $data_dir
fi
#-------------------------


./analysis -n $n -output_dir ${outputdir} -input_path ${inputpath} -ext ${input_ext}

