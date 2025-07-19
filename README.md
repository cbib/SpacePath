[![PyPI - Python Version](https://img.shields.io/pypi/v/smpath)](xxxxx)

smpath
===
A library for spatial metabolomics pathway-based exploration



# Introduction

The `smpath` library facilitates the construction of pathway-based workflows in Python, dedicated to annotated Spatial Metabolomics data.

When working on spatial metabolomics bioinformatic pipelines, **smpath** can be 
easily articulated with other popular Python libraries in the spatial omics world, such as `squidpy`, `scanpy`, `anndata`.
💡

**What is inside?**

Within `smpath`, the computations dedicated to pathways-based analyses rely on powerful libraries such as `sspa`,  `rdkit` and `fastcluster`.

____

_Note_: Currently, the library supports as input data, the aligned spectra datasets and their annotations obtained via [METASPACE](https://metaspace2020.org/). We provide auxiliary scripts for data downloading and formatting, powered by `metaspace-converter` and `pyimzml`.



# Installing smpath 

## System Requirements
`smpath` installation requires an Unix environment with [python 3.10](http://www.python.org/). 
It was tested under Linux and MacOS environments.


## Installation 

The full installation process should take less than 15 minutes on a standard computer.

Via pip command:
`pip install smpath`

One dependency, the `sspa` package, is required to be installed by the user,
 by running `pip install sspa==1.0.4`

All the install process must occur in a virtual environment (e.g.  `
micromamba activate smpath`).

<details><!--comment:micromamba-->

 <summary>Explanation for creating a virtual environment
 <sub><sup>(click to show/hide)</sup></sub>
 </summary>

The users can employ the provided <code>smpath/tools/smpath_env.yaml</code> file:
```
micromamba create -f smpath_env.yaml  
```
or they can alternatively run:
```
micromamba create -n smpath python=3.10.18 -c conda-forge poetry==2.1.3
```

whether using conda or mamba, the arguments after the 'create' command are the same.
Mamba or micromamba are faster, for further information see
 [here](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).
</details>

To sum up, the installation commands in the appropriate order are:

```commandline
micromamba activate smpath  # virtual environment
pip install smpath  # this library
pip install sspa==1.0.4  # required external library 

```

## After installation

For subsequent use just run `micromamba activate smpath`

----

## Developer Setup

To start contributing to smpath, it is required to have a python environment with python >= 3.10 and [poetry](https://python-poetry.org/docs/) installed.

The virtual environment must be activated (see 'Explanation for creating 
of a virtual environment', in the Installation section)

After activating the environment, the project can be installed with
```bash
poetry install
```
Finally install sspa `pip install sspa==1.0.4`

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
