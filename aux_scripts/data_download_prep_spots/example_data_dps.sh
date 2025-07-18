#!/bin/bash

input_file="example-datasets-parameters.tsv"
input_directory="../../datasets_parameters"

data_dir="../../../examples-sm-data"
figures_dir="../../../figures"


echo "Downloading datasets"
python3 metaspace_data_download.py $input_file $input_directory $data_dir

echo ""
echo "Running spot quality classification for all tissues"
python3 spots_purge.py $input_file $data_dir $figures_dir

echo "done"