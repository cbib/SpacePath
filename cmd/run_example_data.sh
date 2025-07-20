#!/bin/bash
# Instruction: must be launched from mymethod/ folder

input_file="datasets_parameters/example-datasets-parameters.tsv"

data_dir="../examples-sm-data"
figures_dir="../figures"


# start
python3 -m mymethod $input_file $data_dir $figures_dir
