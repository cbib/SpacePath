#!/bin/bash
# Instruction: must be launched from mymethod/ folder

input_file="datasets_parameters/public-datasets-parameters-round1.tsv"

data_dir="../data"
figures_dir="../figures"

# start
python3 -m mymethod $input_file $data_dir $figures_dir


