# All what is needed for running spatial metabolomics scripts and notebooks


## 0. Slurm interactive in BB8

Do
```
module load sinteractive
sinteractive --nodes=1 --ntasks-per-node=8
```
Identify the `cortex??` where `??` is the cortex number that was allocated, 
this is the `$MY_NODE` to replace in the step **2.**.





Note: alternatively to Slurm interactive, I can also use a .sh file to launch with sbatch; this is preferred if using a script.


## 1. Environments for spatial metabolomics

**first of all** go to the adequate folder (the parent directory of code repo, data and output):

```
cd /mnt/cbib/thesis_gbm/spatial_thesis/spatial_metabolomics
```

### Python environment
The conda environment  was installed by Aurélien into the module `Metabo/stable`.

It is automatically activated if running in bb8:

```
module load Metabo/stable

```

Note: the .yaml file that Aurélien used is `prototype2_py.yml`


### R environment (if needed)

**Only if GSVA** is necessary (e.g. called in a given R script), run `conda activate r_spatial`.  Its .yaml is found in `/mnt/cbib/thesis_gbm/spatial_thesis/r_spatial.yml` 



## 2. Notebook: tunnel

In the interactive slurm session in BB8:
```
jupyter lab --no-browser --ip=0.0.0.0 --port=8886
```

then from my local
```
ssh -t -t jgalvis@bb8.cbib.u-bordeaux.fr -L 8701:127.0.0.1:8884 ssh $MY_NODE -L 8884:localhost:8886
```
and after twice login, go to `localhost:8701` in my local machine

- - - 

Supplemental information


**The environments' binaries locations**

- Python environment:  `/module/apps/miniconda/python3_isilon/envs/prototype2_py`

- R environment:  /home/jgalvis/.conda/envs/r_spatial

