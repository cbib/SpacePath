 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
 
METHOD
===

METHOD: A pathway-guided framework for deep exploration of Spatial Metabolomics and Lipidomics data



# Introduction

METHOD is a novel bioinformatics framework designed to analyze spatial metabolomics and lipidomics in the context of metabolic pathways. It addresses challenges like compound identification ambiguities and incomplete pathway databases, enabling the mapping of metabolic processes onto tissue architecture.

METHOD's approach integrates pathway activity scores at the spatial resolution level, and investigates unmapped features with the help of molecular similarity approaches. It overcomes limitations of existing tools that lack spatial resolution capabilities, providing a generalizable solution for tissue-specific metabolic mapping, and facilitating the interpretation of complex mass spectrometry imaging data. When applied to cancer tissues, METHOD can enhance the understanding of tumor metabolism.

<img width="900" height="570" alt="Screenshot from 2025-07-22 15-06-39" src="https://github.com/user-attachments/assets/5a5bcc3b-d561-4ffc-b110-5fe0b40b75ca" />

Notably, METHOD was applied on a collection of different types of healthy and tumoral tissues, including glioblastoma. By integrating complex spatially resolved lipidomics data and functional sets expressed as pathways or families, the method provided valuable insights about the glioblastoma's metabolic organization and potential biomarker discovery.
[REF paper or biorxiv XXXXX]()


## How METHOD is implemented

Our method is implemented as a pipeline, that relies on the new `smpath` library. The `smpath` library [link](), offers rich functions for the pathway-based analysis of annotated Spatial Metabolomics data; by implementing `sspa`, `rdkit`, `fastcluster`, among other powerful libraries, for relevant steps.

The METHOD pipeline can also be re-used and adapted to the needs of any type of research team working with Spatial Metabolomics or Spatial Lipidomics data, being interested in the
functional interpretation of their data oriented by the metabolic pathways or families.

_Note_: Currently, our method supports as input data, the aligned spectra datasets and their annotations obtained via [METASPACE](https://metaspace2020.org/). 


# Using METHOD 

## System Requirements
Running METHOD requires an Unix environment with [python 3.10](http://www.python.org/). 
It was tested under Linux and MacOS environments.


## Installation 

The full installation process should take less than 15 minutes on a standard computer. 

First **clone this repository** in your workspace, 

Then install the two dependencies: smpath and sspa. Visit  [LINK TO SMPATH REPO]()
```
pip install smpath
pip install sspa==1.0.4
```
It is strongly recommended to use a virtual environment. [LINK TO SMPATH REPO]()


 
<details>

 <summary>Running unitary tests <sub><sup>(click to show/hide)</sup></sub>
 </summary>
 
We include some unitary tests (this pipeline is a project in active development). 
Run `pytest` from the terminal.

</details>

 
## After installation

For subsequent use just run `micromamba activate smpath` before launching a method command

----

# Documentation

Documentation and **examples of how to use our pipeline** are provided in the [Wiki]() page.   

-----------------------------------------------


## Code organization

I propose to put the code organization into the Wiki

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


