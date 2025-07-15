#!/bin/bash

input_file="example-datasets-parameters.tsv"
input_directory="../../datasets_parameters"

data_dir="../../../data"

# start
python3 metaspace_data_download.py $input_file $input_directory $data_dir