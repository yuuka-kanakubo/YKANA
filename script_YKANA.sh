#!/bin/sh

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


 ## --CentralityCut 2 \
 ## --CentralityCutsumEt \
 ## -n 10000 \
 ## --shuffle \
 ## --BSTR 5 \
 ## --hist_ZeroCentered
./YKANA  \
 -n 10000 \
 -dir /n/work02/yukanaku/mcaa/data/18Apr2024_wohs_K2.0sat2.0 \
 -outdir dndcoordy20_30_18Apr2024_wohs_K2.0sat2.0  \
 -f ev  \
 -ext jets_20_30_pbpb5020.dat  \
 -obs dndcoordy  \
 --EKRTformat  \
 --EKRTbinary \
 --hist_ZeroCentered
