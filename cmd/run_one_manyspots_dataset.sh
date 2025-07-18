#!/bin/bash

input_file="datasets_parameters/manyspots-one-dataset.tsv"

data_dir="../example-manyspots"
figures_dir="../figures"

# start
python3 main.py $input_file $data_dir $figures_dir
