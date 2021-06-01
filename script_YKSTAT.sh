#!/bin/sh

#DIR="../DCCI/data/20210511_PBPB_MB_0Kto5K"
#DIR="../pythia8244/default_pythia_mymain/data/20210501_DEFAULT_PBPB2760GEV_weakStop_100K_sigma0DecayOn"
DIR="../pythia8244/default_pythia_mymain/data/20210426_DEFAULT_PP7TEV_weakStop_1M_sigma0DecayOn"
EV="ev"
#EXT="hadronFinal_corecorona_weakStop.txt"
EXT="default.txt"
outputdir="RtSpectra_PP_PYTHIA"
n=1000


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

./YKSTAT \
 -n $n -outdir ${outputdir} -dir ${DIR} -f ${EV} -ext ${EXT} \
 > ${log_dname}${log_fname}.log  2>&1  &
