#!/bin/sh

DIR="../pythia8244/default_pythia_mymain/data/20210325_DEFAULT_PBPB2760GEV_weakStop_100K"
EV="ev"
#EXT="hadronFinal_core_weakStop_wcol.txt"
EXT="default.txt"
outputdir="PYTHIA_Default_PBPBC24Nch_CORE_WSCAT"
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
# --tagged \
# --2PCfull \
# --2PCnearside \
# --2PCout \
# --setNcoeff 3 \
# --only_corona_associates \
# --vs_Multi 2 \

./YKANA \
 -n $n -outdir ${outputdir} -dir ${DIR} -f ${EV} -ext ${EXT} \
 --4particle \
 --HI \
 > ${log_dname}${log_fname}.log  2>&1  &
