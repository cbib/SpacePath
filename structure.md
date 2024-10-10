# Structure to reproduce the workflow

The workflow uses spatial metabolomics data as input, as well as reference data files, that must be organized exactly as follows.

## required global structure OUTSIDE this repo

```
my_spatial_metabo_project 
    ├── data
    └── output

```
### content of data and output

#### Data
The `data` folder must contain:

* the experimental data: one folder for each acquisition, with a good folder name (in our example below the "brain-ellis" folder is an example). Each acquisition must have its imzML file, ibd file, and metaspace annotations downloaded without modifications.
  
* the reference data: corresponds to two types of files:
  * The gmt file(s) that we need to use. A gmt file contains a list of pathways and their sets of molecules' identifiers within each pathway. In our workflow, these molecules are the metabolites' KEGG identifiers. Several gmt file(s) can be added, making sure that pathways names are unique across all files (should be merged programmatically).
  * The identifiers as they appear in HMDB database, organized in tables, and containing their correspondances to KEGG identifiers. Other tables (LIPIDMAPS, SWISSPROT, etc) can be also saved here. To reproduce our workflow, the idenfitiers table(s) must allways map to KEGG identifiers because the GMT file is given with KEGG identifiers. We provide a table with HMDB identifiers to KEGG ones, and another table with LIPIDMAPS to KEGG identifiers.

  
```
data
├── brain-ellis ---> one MSI acquisition of normal brain (from Metaspace)
│   ├── brain-ellis_hmdb_annotations.csv
│   ├── brain-ellis_lipidmaps_annotations.csv
│   ├── brain-ellis_pos_20um.ibd
│   └── brain-ellis_pos_20um.imzML
│    
└── reference ---> tables containing reference data pulled from databases with external scripts and/or packages
    ├── ..
    ├── ..
    ├── gmt1   <- the KEGG gmt file (generated with ssPA package)
    └── gmt2  <-  a custom gmt file (metabolites in our workflow) given by the biologist expertise

```

#### Output
The `output` folder must contain:

```
output
  tables
  objects
  figures
```
