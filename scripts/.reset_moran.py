import os
import anndata as ad
adata_mz = ad.read_h5ad(filename=os.path.join(
    "../../examples-sm-data", "LP-example0", f'LP-example0.h5ad'))
del adata_mz.uns['moranI']
adata_mz.write_h5ad(filename=os.path.join(
    "../../examples-sm-data", "LP-example0", f'LP-example0.h5ad'))