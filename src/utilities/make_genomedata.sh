#!/bin/bash
# This script gets path of directory of assays wigfix (-a), SNIPER LVs wigfix (-l), genome size file (-g)
 #and create genomedata of those signals

while getopts a:l:g:n: option
do
    case "${option}"
        in
        a) assays_data_path=${OPTARG};;
        l) lv_data_path=${OPTARG};;
        g) chr_size_file_path=${OPTARG};;
	n) file_name=${OPTARG};;
    esac
done

cmd="genomedata-load -s $chr_size_file_path"
#for file in "$assays_data_path"/*.wigfix "$lv_data_path"/*.wigfix ; do
#for file in "$assays_data_path"/*.wigfix; do
for file in "$lv_data_path"/*.wigfix; do
    filename=$(basename "$file")
    if [[ ! $filename = LV* ]]; then
        filename=${filename%.fc*}
        filename=${filename#*-}
    else
        filename="${filename%.*}"
    fi
    cmd="$cmd -t $filename=$file"
done
cmd="$cmd --sizes "$file_name""
echo $cmd
eval "$cmd"
