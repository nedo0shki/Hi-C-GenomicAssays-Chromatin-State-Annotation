#!/bin/bash
# This script gets path of directory of bigwigs (-d), resolution (-r)
 #and create wigfix of those files with a specified resolution

while getopts d:r: option
do
    case "${option}"
        in
        d) DIR_PATH=${OPTARG};;
        r) RES=${OPTARG};;
    esac
done

for path in "$DIR_PATH"/*.bigwig
do
    filepath=`dirname "$path"`
    filename=`basename "$path"`
    filename="${filename%.*}"
    out_filename="${filename}.wig"
    out_filepath="$filepath/$out_filename"
    wigfix_filepath="$filepath/$filename-VirRes$RES.wigfix"
    if test -f "$out_filepath"; then
        echo "$out_filepath exist"
        if test -f "$wigfix_filepath"; then
            echo "$wigfix_filepath exist"
        else
	    python -c "import utils; utils.make_vir_res('$out_filepath', '$RES')"
            rm "$out_filepath"
        fi
    else
        echo $wigfix_filepath
        if test -f "$wigfix_filepath"; then
            echo "$wigfix_filepath exist"
        else
            echo "Making wig file of $filename ..."
            ./bigWigToWig "$path" "$out_filepath"
            echo "Making wigfix file of $filename ..."
	    python -c "import utils; utils.make_vir_res('$out_filepath', '$RES')"
            rm "$out_filepath"
        fi
    fi
done
