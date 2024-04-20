#!/bin/bash

#SBATCH --job-name=ykana
#SBATCH -e errorout
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00

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

./YKANA  -n $1 -outdir $2 -dir $3 -f $4 -ext $5  --EKRTformat  --EKRTbinary > log/$2YKANA.log 
