#!/bin/bash
# Instruction: must be launched from SPACEPATH/ folder

input_file="datasets_parameters/example-datasets-parameters.tsv"

data_dir="../examples-sm-data"
figures_dir="../examples-figures"


# start
python3 -m SPACEPATH $input_file $data_dir $figures_dir
