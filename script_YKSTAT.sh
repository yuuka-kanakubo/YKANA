#!/bin/sh

DIR="../DCCI/data/PP5TEV_transSmear0.6_sigma0.3_pT01.9"
EV="ev"
EXT="hadronFinal_corecorona_K0decay_wcol.txt"
#EXT="default.txt"
outputdir="PP5_MEANPT_K0S"
n=850000


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

./YKSTAT \
 -n $n -outdir ${outputdir} -dir ${DIR} -f ${EV} -ext ${EXT} \
 --vs_dNdeta 10 \
 --ID 310 \
 --CentralityCut_ext hadronFinal_corecorona_weakStop_wcol.txt \
 --INEL_lg_0 \
 > ${log_dname}${log_fname}.log  2>&1  &
