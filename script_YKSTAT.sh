#!/bin/sh

#DIR="../DCCI/data/20210324_PBPB_PT0REF0.9_SIGMA0.5_tau0Unlimited_for_DCCI_CF_T165_MB_0Kto10K_NOparticlization"
#DIR="../DCCI/data/20210322_PP_PT0REF1.8_SIGMA0.5_MB_tau0Unlimited_for_DCCI_CF_T165_0Kto300K_NOparticlization"
DIR="../pythia8244/sample_pythia/data/20210510_JetsTEST"
EV="ev"
EXT="jetinfo.txt"
outputdir="AJ"
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

for i in `seq 0.2 0.2 0.2`
do
for j in `seq 0.5 1.0 0.5`
do
./YKSTAT \
 -n $n -outdir ${outputdir}_a$i"_s"$j -dir ${DIR}_a$i"_s"$j -f ${EV} -ext ${EXT}
# > ${log_dname}${log_fname}_a$i"_s"$j.log  2>&1  &
done
done
