#!/bin/sh

#DIR="../DCCI/data/20210212_PBPB_PT0REF0.9_SIGMA0.5_tau0Unlimited_for_DCCI_CF_T165_NOparticlization"
DIR="../DCCI/data/20210302_PP_PT0REF1.8_SIGMA0.5_tau0Unlimited_for_DCCI_CF_T165_MB_0Kto100K_NOparticlization"
#DIR="../pythia8244/default_pythia_mymain/data/20210305_DEFAULT_PP13TEV"
EV="ev"
EXT="hadronFinal_corecorona_weakStop_wcol.txt"
#EXT="default.txt"
outputdir="4PARTMULTI_PP_DCCI_3sub_ETA0.8"
n=100000


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
 --4particle \
 --threesub  \
 > ${log_dname}${log_fname}.log  2>&1  &

