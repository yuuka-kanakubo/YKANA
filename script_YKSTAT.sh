#!/bin/sh

DIR="../DCCI/data/PBPB2760GEV_transSmear0.6_sigma0.3_pT01.0"
#DIR="../pythia8244/default_pythia_mymain/data/20210305_DEFAULT_PP13TEV"
#DIR="../pythia8244/default_pythia_mymain/data/20210501_DEFAULT_PBPB2760GEV_weakStop_100K_sigma0DecayOn"
#DIR="../pythia8244/default_pythia_mymain/data/20210426_DEFAULT_PP7TEV_weakStop_1M_sigma0DecayOn"
#DIR="../pythia8244/highpt_mode_pythia/data/20210602PP7TEV_HIGHPT_100KxnBin"
EV="ev"
EXT="hadronFinal_corecorona_weakStop.txt"
#EXT="default.txt"
outputdir="V2PT_PBPB2760GEV_transSmear0.6_sigma0.3_pT01.0_105K_WOCOL_CORONA"
n=105000


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
 --only_corona \
 --HI \
 --twosub  \
 --xaxis 3 \
 --xaxis3_input input/BinSettings_V2PT_PBPB.txt \
 --CentralityCut 2 \
 > ${log_dname}${log_fname}.log  2>&1  &
