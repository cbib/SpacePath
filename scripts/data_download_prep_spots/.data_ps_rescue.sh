#!/bin/bash

input_files=("public-datasets-parameters.tsv" "private-datasets-parameters.tsv")
input_directory="../../datasets_parameters"

data_dir="../../../data"
figures_dir="../../../figures"

for input_file in ${input_files[@]}; do
  echo $input_file

  echo "Running spots classification, and if specified, discarding the spots of low quality"
  echo " note: assumes that data_download.sh successfully downloaded the datasets"
  echo ""

  echo "Running custom spot quality classification for TD tissues only"
  echo "  this can be slow, answer 'n' or 'no' if already performed"
  python3 spots_preclassif_TD.py  # no args
  echo ""
  echo "Running spot quality classification for all tissues"
  python3 spots_purge.py $input_file $input_directory $data_dir $figures_dir
  echo "done"
  done