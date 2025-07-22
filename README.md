METHOD
===

METHOD: A pathway-guided framework for deep exploration of Spatial Metabolomics and Lipidomics data



# Introduction

METHOD is a novel bioinformatics framework designed to analyze spatial metabolomics and lipidomics data in glioblastoma. It addresses challenges like compound identification ambiguities and incomplete pathway databases by employing probabilistic methods, TIC normalization, and spatial clustering. This enables accurate mapping of metabolic processes onto tissue architecture, enhancing understanding of tumor metabolism.

METHOD's approach integrates pathway activity scores with spatial clustering, facilitating the interpretation of complex mass spectrometry imaging data. It overcomes limitations of existing tools that lack spatial resolution capabilities, providing a generalizable solution for tissue-specific metabolic mapping. It relies on the new `smpath` library, which offers methods for the pathway-based analysis of annotated Spatial Metabolomics data, available HERE[LINK].

Notably, this method was applied on glioblastoma, and many other types of tissues. Our method offers new insights into glioblastoma's metabolic organization and potential biomarker discovery. 
[REF paper or biorxiv XXXXX]


[DIAGRAM HERE]
💡

## How METHOD is implemented

Our method is implemented as a pipeline, that can be re-used and adapted to the needs of any type of research team working with Spatial Metabolomics or Spatial Lipidomics data and interested in
functional interpretation oriented by metabolic pathways.

_Note_: Currently, our method supports as input data, the aligned spectra datasets and their annotations obtained via [METASPACE](https://metaspace2020.org/). 


# Using METHOD 

## System Requirements
Running METHOD requires an Unix environment with [python 3.10](http://www.python.org/). 
It was tested under Linux and MacOS environments.


## Installation 

The full installation process should take less than 15 minutes on a standard computer. 

First **clone this repository** in your workspace, a recommended folder structure is explained in the Wiki [LINK]

Then install the two dependencies: smpath and sspa. Visit [LINK TO SMPATH REPO]
```
pip install smpath
pip install sspa==1.0.4
```
It is strongly recommended to use a virtual environment  visit [LINK TO SMPATH REPO].

For subsequent uses just run `micromamba activate smpath` before running the method.

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



## After installation

For subsequent use just run `micromamba activate smpath` before launching a method command

----


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


