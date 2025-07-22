 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
 
METHOD
===
A pathway-guided framework for deep exploration of Spatial Metabolomics and Lipidomics data


This repository contains the code for the manuscript "XXXXXXXXXX" [REF paper or biorxiv XXXXX](). The code is also applicable to datasets other than those included in the paper.

# Introduction

METHOD is a novel bioinformatics framework designed to analyze spatial metabolomics and lipidomics in the context of metabolic pathways. It addresses challenges like compound identification ambiguities and incomplete pathway databases, enabling the mapping of metabolic processes onto tissue architecture.

METHOD's approach integrates pathway activity scores at the spatial resolution level, and investigates unmapped features with a strategy that exploits the molecular similarity between chemical entities. It overcomes limitations of existing tools that lack spatial resolution capabilities, providing a generalizable solution for tissue-specific metabolic mapping, and facilitating the interpretation of complex mass spectrometry imaging data. When applied to cancer tissues, METHOD can enhance the understanding of tumor metabolism.

<img width="900" height="570" alt="Screenshot from 2025-07-22 15-06-39" src="https://github.com/user-attachments/assets/5a5bcc3b-d561-4ffc-b110-5fe0b40b75ca" />


Notably, METHOD was applied on a collection of different types of healthy and tumoral tissues, including glioblastoma. By integrating complex spatially resolved lipidomics data and functional sets expressed as pathways or families, the method provided valuable insights about the glioblastoma's metabolic organization and potential biomarker discovery.
[REF paper or biorxiv XXXXX]()

_Note_: Currently, our method supports as input data, the datasets consisting of aligned and recalibrated spectra intensities alongside their annotations, that are all obtained from the [METASPACE](https://metaspace2020.org/) output. 


# Using METHOD 

## System Requirements
Running METHOD requires an Unix environment with [python 3.10](http://www.python.org/). 
It was tested under Linux and MacOS environments.



**_Information about the method implementation_**

Our method is implemented as a pipeline that entirely depends on the new `smpath` library ([link]()), which is very easy to install (see Installation section of this README).

> [!NOTE]  
> Our library `smpath` offers rich functions for the pathway-based analysis of annotated spatial metabolomics data. It implements `sspa`, `rdkit`, `fastcluster`, among other powerful libraries.


## Installation 

The full installation process should take less than 15 minutes on a standard computer (dependencies included). 

1. First **clone `mymethod`** in your workspace (`git clone git@github.com:cbib/mymethod.git`)   %%caution, link will change when name appropriate%%
2. also **clone `smpath`** in your workspace (`git clone git@github.com:cbib/smpath.git`)
3. activate your virtual environment with poetry installed (visit [LINK TO SMPATH REPO]()); 
4. run:
```
cd smpath
poetry install
pip install sspa==1.0.4
```
 
<details>

 <summary>Running unitary tests <sub><sup>(click to show/hide)</sup></sub>
 </summary>
 
We include some unitary tests (this project in active development). 
```
cd mymethod
pytest
```

</details>

 
## After installation

For subsequent use just activate the environment (e.g. `micromamba activate smpath`) as it has all the requirements already installed for running METHOD.
Then you can run the generic command of METHOD (see Use case).


## Use case

The generic command line is:
```
python -m mymethod $INPUT_FILE $DATA_FOLDER $FIGURES_FOLDER
```
----

# Documentation

Documentation, examples of how to use our pipeline, and the specifics of the datasets used in the analysis presented in our manuscript, are provided in the [Wiki]() page.   

-----------------------------------------------


## Code organization


* `src/..  TO COMPLETE HERE !!
* `src/..  TO COMPLETE HERE !!
* `src/..  TO COMPLETE HERE !!
* `tests` directory contains unit tests


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


