#!/bin/bash


input_files_list=("public-datasets-parameters.tsv" "private-datasets-parameters.tsv")
input_directory="../../datasets_parameters"

data_dir="../../../data"
figures_dir="../../../figures"

for i_file in ${input_files_list[@]}; do

  input_file=${input_directory}/${i_file}

  echo $input_file

  echo "Downloading datasets"
  python3 metaspace_data_download.py $input_file $input_directory $data_dir

  echo ""

  echo "Running custom spot quality classification, for selected tissues only"
  echo "  (TD tissues) this can be slow, answer 'n' or 'no' if already performed"
  python3 spots_preclassif_TD.py  # no args

  echo " (brain 7 tissue) if already performed, answer 'n'"
  python3 spots_preclassif_LP7.py


  echo ""
  echo "Running spots classification, and if specified, discarding the spots of low quality"

  echo "Running spot quality classification for all tissues"
  python3 spots_purge.py $input_file $data_dir $figures_dir
  echo "done"

  done