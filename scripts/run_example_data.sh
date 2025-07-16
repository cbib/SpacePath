#!/bin/bash

input_file="example-datasets-parameters.tsv"
input_directory="../datasets_parameters"

data_dir="../../examples-sm-data"
figures_dir="../../figures"


# start
python3 main.py $input_file $input_directory $data_dir $figures_dir
