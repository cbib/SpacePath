# smpath

A comprehensive library to construct a pathway-based analysis of annotated Spatial Metabolomics data.

**smpath** powers the method XXXXX, available at @githubXXXXX. 
💡

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

Developers can easily install a local version with poetry.
Clone the repository and then run:
```
cd smpath
poetry install
```
</details>

# License


TO FILL
