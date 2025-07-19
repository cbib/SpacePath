# smpath

A comprehensive library to construct a pathway-based analysis of annotated Spatial Metabolomics data.


In the context of spatial metabolomics bioinformatic pipelines, **smpath** can be easily articulated with other  popular libraries in the spatial omics world, such as `squidpy`, `scanpy`, `anndata`.
💡

**What is inside?**

Within `smpath`, the modules dedicated to pathways-based analyses rely on powerful libraries such as `sspa`,  `rdkit` and `fastcluster`.


_Note_: The aligned spectra datasets and their annotations obtained via [METASPACE](https://metaspace2020.org/) are currently supported. We also provide auxiliary scripts for data downloading and formatting, powered by `metaspace-converter` and `pyimzml`.

# Requirements

smpath runs under Python >= 3.10. It is strongly recommended to install in a mamba, micromamba or conda virtual environment
(see [here](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)).



<details><!--comment:micromamba-->

 <summary>Explanation of the creation of a micromamba environment
 <sub><sup>(click to show/hide)</sup></sub>
 </summary>

The users can employ a .yaml file, as the one provided in <code>smpath/tools</code>
```
micromamba create -f mambaenv.yaml  
```
or they can alternatively run:
```
micromamba create -n mambaenv python=3.10.18 -c conda-forge poetry==2.1.3
```
If the user has a preference for conda or miniconda, the word micromamba must be replaced by conda or miniconda, but the rest of the line after the 'create' command is the same.
</details>

Once the virtual environment is created, activate it
```
micromamba activate mambaenv
```

# Install

Users can easily install `smpath` via pip

```
pip install smpath
```

<details><!--comment:micromamba-->

 <summary>For developers
 <sub><sup>(click to show/hide)</sup></sub>
 </summary>

After setting up and activating the virtual environment, clone the repository, locate yourself inside the cloned directory, and then run poetry:
```
micromamba activate mambaenv
cd smpath
poetry install
```
</details>

# Documentation

Documentation and examples of how to use with `smpath` library are provided in the Wiki page.

# License


TO FILL
