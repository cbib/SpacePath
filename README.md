 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
 
SPACEPATH
===
A pathway-guided framework for deep exploration of Spatial Metabolomics and Lipidomics data


This repository contains the code of the method SPACEPATH (Spatial omics exploration based on molecular families and metabolic pathways profiles), dedicated to the analysis of annotated Spatial Metabolomics and Lipidomics data. The code is also applicable to annotated datasets other than those included in our study.

# Introduction

SPACEPATH is a novel bioinformatics framework designed to analyze spatial metabolomics and lipidomics in the context of metabolic pathways. It addresses challenges like compound identification ambiguities and incomplete pathway databases, enabling the mapping of metabolic processes onto tissue architecture.

SPACEPATH's approach integrates pathway activity scores at the spatial resolution level, and investigates unmapped features with a strategy that exploits the molecular similarity between chemical entities. It overcomes limitations of existing tools that lack spatial resolution capabilities, providing a generalizable solution for tissue-specific metabolic mapping, and facilitating the interpretation of complex mass spectrometry imaging data. When applied to cancer tissues, SPACEPATH can enhance the understanding of tumor metabolism.

<img width="900" height="570" alt="Screenshot from 2025-07-22 15-06-39" src="https://github.com/user-attachments/assets/5a5bcc3b-d561-4ffc-b110-5fe0b40b75ca" />


Notably, SPACEPATH was applied on a collection of different types of healthy and tumoral tissues, including glioblastoma. By integrating complex spatially resolved lipidomics data and functional sets expressed as pathways or families, the method provided valuable insights about the glioblastoma's metabolic organization, prompting a potential impact on biomarker discovery.

_Note_: Currently, our method supports as input data, the datasets consisting of aligned and recalibrated spectra intensities alongside their annotations, that are obtained from the [METASPACE](https://metaspace2020.org/) output. 


# Installing SPACEPATH 

## System Requirements
Running SPACEPATH requires an Unix environment with [python 3.10](http://www.python.org/). 
It was tested under Linux and MacOS environments.

**_Information about the method implementation_**: SPACEPATH is a Pipeline that entirely depends on the new **`smpath`** library, which is very easy to install (see how to [Install](https://github.com/cbib/smpath?tab=readme-ov-file#installing-smpath) this core library).

> [!NOTE]  
> The associated core library `smpath` offers rich functions for the pathway-based analysis of annotated spatial metabolomics data, and for exploiting unmapped features. It implements `sspa`, `rdkit`, `fastcluster`, among other powerful libraries.


## Installing SPACEPATH 

The full installation process should take less than 15 minutes on a standard computer (dependencies included). 

1. First **clone `SPACEPATH`** in your workspace (`git clone git@github.com:cbib/SPACEPATH.git`) 
2. also **clone `smpath`** in your workspace (`git clone git@github.com:cbib/smpath.git`)
3. activate your virtual environment with poetry installed (visit [link](https://github.com/cbib/smpath)); 
4. run:
```
cd smpath
poetry install
pip install sspa==1.0.4
```
_Note_: the associated core library `smpath` which will be soon be packaged into PyPI to be installable via `pip`. 
<details>

 <summary>Running unitary tests <sub><sup>(click to show/hide)</sup></sub>
 </summary>
 
We include some unitary tests in this repository, but in essence this code is a Pipeline and not a package. The most relevant unitary tests are implemented in the core library `smpath`. 
```
cd SMPATH
pytest
```
</details>
 
## After the install process

For subsequent use, just activate the environment (e.g. `micromamba activate smpath`) as it has all the requirements already installed for running SPACEPATH.

# Using SMPATH

The generic command line is:
```
python -m mymethod $INPUT_FILE $DATA_FOLDER $FIGURES_FOLDER
```
----

# Documentation

Documentation, examples of how to use our pipeline, and the specifics of the datasets used in the analyses 
are provided in the [Wiki]() page.   

-----------------------------------------------------------------------------------------------
  
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


