#!/bin/sh


inputpath="../DCCI/data/20201221_PBPB_PT0REF0.9_SIGMA0.5_PTCUT_minusINF"
input_ext="hadronFinal_corecorona_weakStop.txt"
outputdir="VERTICES_PBPB_PT0REF0.9_SIGMA0.5_PTCUT_minusINF"
n=10

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


./analysis -n $n -outdir ${outputdir} -path ${inputpath} -ext ${input_ext}

