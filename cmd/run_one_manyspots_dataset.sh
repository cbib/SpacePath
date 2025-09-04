#!/bin/bash
# Instruction: must be launched from SPACEPATH/ folder

input_file="datasets_parameters/manyspots-one-dataset.tsv"

data_dir="../example-manyspots"
figures_dir="../figures"

# start
python3 -m SPACEPATH $input_file $data_dir $figures_dir
