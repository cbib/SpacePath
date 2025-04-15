**all the following is no longer valid because the "whole-ellis" data 
was qTOF and so no suitable for pathways analyses**



1. Note that for the "whole-ellis" data and its scripts, that the **aux** "zonist" script: `matrix_from_metaspace.py`
which is at https://github.com/johaGL/zonist/trese/master/src/zonist/aux/
must be used for generating the `whole-ellis_mx.h5ad` file !!!!  --> the new script to get-adata was done, and the scripts were renamed. OK

2. the figures are at `user-pipeline/scripts/whole-ellis/0-figures` but it is not a good location !!! --->  relocated to figures/, OK

3. the notebooks , specifically in the sections of the data MATRIX retrieval,  must explain that: a)  "brain-ellis" is a small subsample for giving a working example. b) ways of getting the MATRIX are three: (i) download  from METASPACE website (button), (ii) from imzML+ibd or (iii) use metaspace_converter library.  First organize the "zonist" scripts that do (ii) and (iii), and after that, give reference to them in the notebooks 
