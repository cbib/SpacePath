# smpath

A comprehensive library to construct a pathway-based analysis of annotated Spatial Metabolomics data.



In the context of spatial metabolomics bioinformatic pipelines, **smpath** can be easily articulated with other  popular libraries in the spatial omics world, such as `squidpy`, `scanpy`, `anndata`.
💡

**What is inside?**

Within `smpath`, the modules dedicated to pathways-based analyses rely on powerful libraries such as `sspa`,  `rdkit` and `fastcluster`.


_Note_: The aligned spectra datasets and their annotations obtained via [METASPACE]() are currently supported. We also provide auxiliary scripts for data downloading and formatting, powered by `metaspace-converter` and `pyimzml`.

# Requirements

smpath runs under Python >= 3.10. It is strongly recommended to install in a mamba, micromamba or conda virtual environment
see [here](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).

%%%%Note for myself: do as here: https://blog.stackademic.com/conda-and-poetry-a-harmonious-fusion-8116895b6380 %%%%

<details><!--comment:micromamba-->

 <summary>Explanation of the creation of a micromamba environment
 <sub><sup>(click to show/hide)</sup></sub>
 </summary>

Users can employ a .yaml file, as the one provided in <code>smpath/tools</code>
```
micromamba create -f mambaenv.yml  
```
or they can alternatively create it with:
```
micromamba create -n mambaenv python=3.10 -c conda-forge
```
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
