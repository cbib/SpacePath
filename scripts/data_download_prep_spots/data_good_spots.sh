#!/bin/bash

input_file="example-datasets-parameters.tsv"  # "private-datasets-parameters.tsv"
input_directory="../../datasets_parameters"

data_dir="../../../data"
figures_dir="../../../figures"

echo "Preparing dat at the level of spots by classifying spots by their quality, and if specified, discarding the spots of low quality"
echo " note: if not done, please run data_download.sh in first place"
echo ""

echo "Running custom spot quality classification for TD tissues only"
echo "  this can be slow, answer 'n' or 'no' if already performed"
python3 spots_preclassif_TD.py  # no args
echo ""
echo "Running spot quality classification for all tissues"
python3 spots_purge_final.py $input_file $input_directory $data_dir $figures_dir