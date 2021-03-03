#!/bin/sh


DIR="../DCCI/data/20210129_PP_PT0REF1.8_SIGMA0.5_tau0Unlimited_for_DCCI_CF_T165_0Kto50K_NOparticlization"
EV="ev"
EXT="hadronFinal_corecorona_weakStop.txt"
outputdir="PP7TEV_DCCI_MEANMT_WOCOL"
n=50000


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


./analysis -n $n -outdir ${outputdir} -dir ${DIR} -f ${EV} -ext ${EXT} --INEL_lg_0 --CentralityCut 9 > ${log_dname}${log_fname}.log   &

