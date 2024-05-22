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

##./YKANA  -n $1 -outdir $2 -dir $3 -f $4 -ext $5  \
## --CentralityCut 2 \
## -print_dndeta \
## --xaxis 1 \
## --EKRTformat  \
## --EKRTbinary > log/$2YKANA.log 

./YKANA  \
 -n 1000 \
 -dir /n/work02/yukanaku/mcaa/data/30Apr2024_whs_K2.4sat2.0_allev_1000  \
 -outdir CentralityCutswhs_K2.4sat2.0_1000  \
 -f ev  \
 -ext jets_0_100_pbpb5020.dat  \
 --CentralityCut 2 \
 --shuffle \
 --EKRTformat  \
 --EKRTbinary \
 --xaxis 1 
