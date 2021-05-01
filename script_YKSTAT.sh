#!/bin/sh

DIR="../DCCI/data/DCCI_PBPB_MB_0to10K"
EV="ev"
EXT="hadronFinal_corecorona_weakStop.txt"
outputdir="MTSCALING_PBPB_DCCI_CORONA_WOCOL"
n=10000


#Do not modify this.
#---------------------
log_dname="log/"
data_dir="data/"

today=$(date "+%Y%m%d")
log_fname=${today}${outputdir}
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


#Options.
#--------------
# --CentralityCut 9 \
# --CentralityCut_ext hadronFinal_corecorona_weakStop.txt \
# --HI \
# --only_corona \
# --INEL_lg_0 \
# --twosub  \
# --threesub  \
# --4particle \

./YKSTAT \
 -n $n -outdir ${outputdir} -dir ${DIR} -f ${EV} -ext ${EXT} \
 --HI \
 --only_corona \
 --CentralityCut 9 \
 > ${log_dname}${log_fname}.log  2>&1  &
