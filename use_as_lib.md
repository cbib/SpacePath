# How to use "zonist" as a library

In Server:


Locate yourself in the parent folder of both data folder and zonist local clone. In our case the parent folder is "spatial_metabolomics/"

* spatial_metabolomics
    * data
    * output
    * zonist  (cloned repository)


```
cd spatial_metabolomics
```


Activate an interactive Slurm session
```
module load sinteractive
sinteractive --nodes=1 --ntasks-per-node=8

```

Activate the conda environment, launch python:
```
module load Metabo/stable
python
```


```
from zonist.src import spectra2anndata

spectra = spectra2anndata.open_spectra("my_imzml","my_ibd")
```

and call all the rest of functions you want


finally, close slurm connection
```
scancel <jobid>
```

