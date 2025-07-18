import os
import anndata as ad
SAMPLE = "LP-example0"
adata_mz = ad.read_h5ad(filename=os.path.join(
    "../../examples-sm-data", SAMPLE, f'{SAMPLE}.h5ad'))
del adata_mz.uns['moranI']
adata_mz.write_h5ad(filename=os.path.join(
    "../../examples-sm-data", SAMPLE, f'{SAMPLE}.h5ad'))