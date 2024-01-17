#!/bin/sh

#SBATCH --job-name=allon_sum
#SBATCH -e errorout
#SBATCH --output=log/outYKANA.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00


DIR="/n/work02/yukanaku/mcaa-master/data/16Jan2024_test/60_70"
EV="ev"
##EXT="hadronFinal_core_weakStop_wcol.txt"
EXT="EKRTminijet.txt"
outputdir="test"
n=1


##Do not modify this.
##---------------------
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
##-------------------------


##Options.
##--------------
## --CentralityCut 9 \
## --CentralityCut_ext hadronFinal_corecorona_weakStop.txt \
## --HI \
## --only_corona \
## --INEL_lg_0 \
## --twosub  \
## --threesub  \
## --4particle \
## --tagged \
## --2PCfull \
## --2PCnearside \
## --2PCout \
## --setNcoeff 3 \
## --only_corona_associates \
## --vs_Multi 2 \

./YKANA \
 -n $n -outdir ${outputdir} -dir ${DIR} -f ${EV} -ext ${EXT} \
 --EKRTformat 
 ##> ${log_dname}${log_fname}.log  2>&1  &
