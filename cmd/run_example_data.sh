#!/bin/bash
# Instruction: must be launched from SpacePath/ folder

input_file="datasets_parameters/example-datasets-parameters.tsv"

data_dir="../examples-sm-data"
figures_dir="../examples-figures"


# start
SpacePath $input_file $data_dir $figures_dir
