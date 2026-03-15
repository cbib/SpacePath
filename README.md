[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


smpath
===
A library that supports the pathway-based and molecular families exploration of annotated spatial metabolomics/lipidomics data

# Introduction

The `smpath` library offers methods for the comprehensive functional-oriented analysis of annotated spatial metabolomics data.
This is the core library of the pipeline [SpacePath](https://github.com/cbib/SpacePath/tree/main)

Moreover, when building any spatial metabolomics bioinformatic pipeline, **smpath** can be 
easily articulated with other popular Python libraries in the _spatial omics_ world, such as `squidpy`, `scanpy`, `anndata` and more 
💡

**What is inside?**

Within `smpath`, the computations dedicated to pathways-based analyses rely on powerful libraries such as `sspa`,  `rdkit` and `fastcluster`.

____

_Note_: Currently, the library supports as input data, the aligned spectra datasets and their annotations obtained via [METASPACE](https://metaspace2020.org/). Auxiliary scripts for data downloading and formatting, powered by `metaspace-converter` and `pyimzml`, are also provided.


# Installing smpath 

## System Requirements
`smpath` installation requires an Unix environment with [python 3.10](http://www.python.org/) . 
It was tested under Linux and MacOS environments.


## Installation 

The full installation process should take less than 15 minutes on a standard computer. A virtual environment with poetry 2.1.3 is required
(see 'Explanation for creating a virtual environment').

<details><!--commentvirtualenv-->

 <summary>Explanation for creating a virtual environment
 <sub><sup>(click to show/hide)</sup></sub>
 </summary>

All the install process must occur in a virtual environment with poetry installed.
* with venv in python (>=3.10,<3.13):
```
python3 -m venv smpath
source smpath/bin/activate  # activate
pip install poetry==2.1.3
```
* alternatiely, with conda, miniconda, micromamba or mamba (command examples below use micromamba):
  
```
micromamba create -n smpath python=3.10.18 -c conda-forge poetry==2.1.3
```
* or use the provided .yaml file that is located at <code>smpath/tools/smpath_env.yaml</code> :
```
micromamba create -f smpath_env.yaml  
```

</details>

Currently `smpath` is installable from this GitHub `git clone git@github.com:cbib/smpath.git`.
The installation commands in the appropriate order are:

```commandline
conda activate smpath  # or `source smpath/bin/activate` if venv, to activate the environment
cd smpath
poetry install   # installs this library
pip install sspa==1.0.4    # required external library 
```


## After installation

For subsequent use just activate your environment (e.g. `micromamba activate smpath`)

----

## Developer Setup

To start contributing to smpath, it is required to have a python environment with python >= 3.10 (and < 3.13) and [poetry](https://python-poetry.org/docs/) installed.

The virtual environment with poetry must be activated (see 'Explanation for creating 
of a virtual environment', in the Installation section)

After activating the environment, run in the project folder:

```bash
poetry install
pip install sspa==1.0.4
```

------

## Code organization

* `src/smpath/processing` directory contains the implemented high-level  functions for the processed dedicated to the measurements and/or the pathways inputs
* `src/smpath/visualization` directory contains the implemented high-level scripts to produce figures
* `src/smpath/data` directory contains the python classes for data initialization alongside the methods proper to each data class
* `src/smpath/method` directory contains the python classes of the configurations and methods that are applied on diverse types of input objects, by invoking 'src/smpath/processing' functions
* `src/smpath/anndata_aux` directory contains the auxiliary functions for dataset download and formatting into an 'anndata' object
* `tests` directory contains unit tests
* `tools` directory contains virtual environment setup .yaml files.


## Running unit tests 

* With pytest, by running `pytest` from `smpath`
* Alternatively, place yourself in `smpath/tests` and execute `python -m unittest` 
* If the project was installed with `poetry install`, tests can also be run using `poetry run pytest` or from VSCode's GUI

-----------------------------------------------------------------------------------------------

# Documentation

Documentation and examples of how to use with `smpath` library 
are provided in the Wiki page.

-----------------------------------------------
  
# Getting help

For any information or help, you can get in touch with: 

* [Johanna Galvis](mailto:deisy-johanna.galvis-rodriguez[AT]u-bordeaux.fr)


# LICENSE MIT

    Copyright (c) 2025
    
    Johanna Galvis (1,2)    deisy-johanna.galvis-rodriguez@u-bordeaux.fr

    Macha Nikolski (1,2)    macha.nikolski@u-bordeaux.fr
    
    (1) CNRS, IBGC - University of Bordeaux,
    1, rue Camille Saint-Saens, Bordeaux, France

    (2) CBiB - University of Bordeaux,
    146, rue Leo Saignat, Bordeaux, France

    
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
