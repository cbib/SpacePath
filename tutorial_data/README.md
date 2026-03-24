# Datasets for using the `smpath` library

Datasets are provided, as examples of spatial metabolomics/lipidomics acquisitions, 
through the public Zenodo registry [https://zenodo.org/records/19113236](https://zenodo.org/records/19113236). 

<div class="alert alert-block alert-info">
<b>Important:</b> . 
These datasets are not contained within the library distribution, instead, they are provided separately as external data for demonstration purposes.
</div>

**Instructions:**

Always locate yourself in a working directory (e.g. `$HOME`):


1. Download the data from 
[https://zenodo.org/records/19113236](https://zenodo.org/records/19113236) 
and unzip the folder.   This is the expected folder structure:

```
$working_dir   # e.g.$HOME or other
├── examples-sm-data      # <-- unzipped data folder
│    ├── example-corn-root      # the provided example dataset N 1
│    ├── metabolites-mappings      # the provided reference pathways
│    ├── example-cerebellum      # the provided example dataset N 2
│
├── smpath    # if local cloning of the library; DO NOT modify
│    
├── my_scripts/    # HERE your .py scripts or jupyter notebooks
```
Details about the specific content of each subfolder of `examples-sm-data` are provided in the same Zenodo registry.


2. Visit our Wiki page [LINK], which shows the usage of `smpath` on 
the provided example data. There we provide code snippets, that you can re-use to run a fast analysis on example data!


__Note__: if you cloned the `smpath` library locally for usage (and not development), 
**do not add any datasets**, nor scripts, inside the library folder.



