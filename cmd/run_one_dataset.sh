#!/bin/bash
# Instruction: must be launched from mymethod/ folder

input_file="datasets_parameters/private-one-dataset.tsv"

data_dir="../data"
figures_dir="../figures"

# start
python3 -m mymethod $input_file $data_dir $figures_dir
