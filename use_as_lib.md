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
cd zonist
pip install -e .


```

Finally do a script in python wherever you want
```

#import sys
# caution: path[0] is reserved for script path (or '' in REPL)
#sys.path.insert(1, '/mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics/zonist/')
#import zonist as zo
#from zo.data import DataSpatialMatrix


# OR
# locate above zonist, launch pip install -e zonist

import zonist.data as zd   # works if pip -e ... thing

##
# vars

mappingfiles = {'HMDB': 'ww', 'LIPIDMAPS': 'mm'} # kegg to kegg is senseless

mx = {'HMDB': 'r', 
      'LIPIDMAPS': 't', 
      'KEGG': ""}  # if no kegg, empty string

annots = {'HMDB': 'oo', 
         'LIPIDMAPS': 'ee',
         'KEGG': ""}  # same as above

##

anno_obj = zd.References(
    parent_dir="nose",    
    subfolder_mappings="mappingsfolder", 
    mappings_filenames=mappingfiles, 
   )
print("anno_obj", anno_obj, "\n")

oh = zd.DataSpatialMatrices(
    name="heheh",
    parent_dir="nose",    
    subfolder_data="heheh",  annotations_filenames=annots,
    matrices_filenames=mx)

print("datatata", oh, "\n")

print("loaded")


```

and call all the rest of functions you want


finally, close slurm connection
```
scancel <jobid>
```

