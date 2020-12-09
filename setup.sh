#!/bin/sh

log_dname="log/"
data_dir="data/"
today=$(date "+%Y%m%d")
log_fname=${today}${dirname}
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
