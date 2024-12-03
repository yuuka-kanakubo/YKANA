#!/bin/sh

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

 ## --CentralityCut 2 \
 ## --sortsumEt \
 ## -n 10000 \
 ## --shuffle \
 ## --BSTR 5 \
 ## --hist_ZeroCentered
 ##--EKRTformat  \
 ##--EKRTbinary \
./YKANA  \
 -n 5000 \
 -dir /home/yukanaku/MC-EKRT/data/20241119_wohsK2.0kappa2.0_pPb/ \
 -outdir MCEKRTpPb5020_0_5_K2.0kappa2.0_detdy_cm_5K \
 -f ev  \
 --EKRTformat \
 --EKRTbinary \
 --pPb_cm2lab \
 --parton \
 -ext jets_0_5_ppb5020.dat  \
 -obs detdy  \
 --hist_ZeroCentered


##CENT=(0_5 20_30 60_70)
##K=(2.0 3.0)
##Kapsat=(1.0 2.0)
##
##for i in `seq 0 1 2`
##do
##	for k in `seq 0 1 0`
##	do
##		for sat in `seq 0 1 1`
##		do
##   		DIR=/n/work02/yukanaku/mcaa/data/18Apr2024_wohs_K${K[k]}sat${Kapsat[sat]}
##   		EV="ev"
##   		##EXT="hadronFinal_core_weakStop_wcol.txt"
##   		EXT=jets_${CENT[i]}_pbpb5020.dat
##   		outputdir=detdy_wohs_K${K[k]}sat${Kapsat[sat]}_${CENT[i]}
##   		n=1000
##   		
##   		
##   		##Do not modify this.
##   		##---------------------
##   		log_dname="log/"
##   		data_dir="data/"
##   		
##   		today=$(date "+%Y%m%d")
##   		log_fname=${today}${outputdir}
##   		if [ ! -d $log_dname ]
##   		then
##   		    echo "Directory "$log_dname" DOES NOT exists." 
##   		    echo "mkdir "$log_dname
##   		    mkdir $log_dname
##   		fi
##   		if [ ! -d $data_dir ]
##   		then
##   		    echo "Directory "$data_dir" DOES NOT exists." 
##   		    mkdir $data_dir
##   		fi
##   		##-------------------------
##   		
##   		sbatch -o YKANA${outputdir}.out  runYKANA.sbatch $n ${outputdir} ${DIR} ${EV} ${EXT}
##   		echo sbatch -o YKANA${outputdir}.out  runYKANA.sbatch $n ${outputdir} ${DIR} ${EV} ${EXT}
##
##		sleep 0.1
##		done
##	sleep 0.1
##	done
##sleep 0.1
##done
##

