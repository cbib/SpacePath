#!/bin/bash
# Instruction: must be launched from SpacePath/ folder


input_file="datasets_parameters/inhouse-datasets-parameters.tsv"

data_dir="../data"
figures_dir="../figures"

# start
SpacePath $input_file $data_dir $figures_dir
