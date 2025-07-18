#!/bin/bash

input_file="datasets_parameters/private-datasets-parameters.tsv"

data_dir="data"
figures_dir="figures"

# start
python3 main.py $input_file $data_dir $figures_dir
