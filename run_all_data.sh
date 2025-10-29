#!/bin/bash
# Instruction: must be launched from SPACEPATH/ folder


input_file="datasets_parameters/private-datasets-parameters.tsv"

data_dir="../data"
figures_dir="../figures"

# start
python3 -m SPACEPATH $input_file $data_dir $figures_dir

input_file2="datasets_parameters/public-datasets-parameters.tsv"

# start
python3 -m SPACEPATH $input_file2 $data_dir $figures_dir

