[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)  [![Poetry](https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
 

SpacePath
===
A framework for pathways and molecular sets exploration of annotated Spatial Metabolomics and Lipidomics data


# Introduction

SpacePath is a novel bioinformatics framework designed to analyze spatial metabolomics and lipidomics in the context of metabolic pathways. It addresses challenges like compound identification ambiguities and incomplete pathway databases, enabling the mapping of metabolic processes onto tissue architecture.

SpacePath's approach integrates pathway activity scores at the spatial resolution level, and investigates unmapped features with a strategy that exploits the molecular similarity between chemical entities. It overcomes limitations of existing tools that lack spatial resolution capabilities, providing a generalizable solution for tissue-specific metabolic mapping, and facilitating the interpretation of complex mass spectrometry imaging data. When applied to cancer tissues, SpacePath can enhance the understanding of tumor metabolism.

<img width="900" height="500" alt="Screenshot from 2025-07-22 15-06-39" src="https://github.com/user-attachments/assets/5a5bcc3b-d561-4ffc-b110-5fe0b40b75ca" />


Notably, SpacePath was applied on a collection of different types of healthy and tumoral tissues, including glioblastoma. By integrating complex spatially resolved lipidomics data and functional sets expressed as pathways or families, the method provided valuable insights about the glioblastoma's metabolic organization, prompting a potential impact on biomarker discovery.

_Note_: Our method also supports, as input data, the datasets consisting of aligned and recalibrated spectra intensities alongside their annotations, that are obtained from the [METASPACE](https://metaspace2020.org/) output. 


# Installing SpacePath 

## System Requirements
Running SpacePath requires an Unix environment with [python 3.10](http://www.python.org/). 
It was tested under Linux and MacOS environments.

 A virtual environment with poetry 2.3.3 is required
(see 'Explanation for creating a virtual environment').

<details><!--commentvirtualenv-->

 <summary>Explanation for creating a virtual environment
 <sub><sup>(click to show/hide)</sup></sub>
 </summary>

All the install process must occur in a virtual environment with poetry installed.
* with venv in python:
```
python3 -m venv spacepath
source spacepath/bin/activate  # activate
pip install poetry==2.3.3
```
* alternatiely, with conda, miniconda, micromamba or mamba (command examples below use micromamba):
  
```
micromamba create -n spacepath python=3.10.18 -c conda-forge poetry==2.3.3
```
* or use the provided .yaml file that is located at <code>tools/spacepath_env.yaml</code> :
```
micromamba create -f spacepath_env.yaml  
```

</details>


## Installation 

The full installation process should take less than 15 minutes on a standard computer (dependencies included). Two options are possible.

#### Option 1: Install via `pip`
Follow the three steps specified in the snippet below:
```
micromamba activate spacepath   # activate your environment
pip install SpacePath         
pip install sspa==1.0.4
```

#### Option 2: Using a local clone and poetry

First clone `SpacePath` in your workspace (`git clone git@github.com:cbib/SpacePath.git`), then, follow the steps specified in the snippet below:
```
micromamba activate spacepath   # 1. activate your environment
cd SpacePath   # 2. locate into the cloned SpacePath
poetry install  # 3. install SpacePath
pip install sspa==1.0.4  # 4. intall sspa
```
 
## After the install process

For subsequent use, just activate the environment (e.g. `micromamba activate spacepath`) as it has all the requirements already installed for running SpacePath.

# Using SpacePath

The generic command line to run the pipeline is:
```
SpacePath $input_file $data_folder $figures_folder
```

# Documentation

Documentation, examples of how to use `SpacePath`, and the specifics of the datasets used in the analyses 
are provided in the [Wiki]() page.   

-------------

# For Developers

### Code organization of the `core` folder of SpacePath

The `core` folder contains the core classes for data and methods that are required for running SpacePath. It is organized into:

* `processing` directory contains the implemented high-level functions for the processed dedicated to the measurements and/or the pathways inputs
* `visualization` directory contains the implemented high-level scripts to produce figures
* `data` directory contains the python classes for data initialization alongside the methods proper to each data class
* `method` directory contains the python classes of the configurations and methods that are applied on diverse types of input objects, by invoking 'core/processing' functions
* `anndata_aux` directory contains the auxiliary functions for dataset download and formatting into an 'anndata' object


### Running unit tests 

Place yourself in `SpacePath/tests` and execute **`pytest`** or,  alternatively, `python -m unittest`. 


-------------
  
# Getting help

For any information or help, you can get in touch with: 

* [Johanna Galvis](mailto:deisy-johanna.galvis-rodriguez[AT]u-bordeaux.fr)


-------------

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


