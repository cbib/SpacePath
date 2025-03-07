import os

import zonist.aux.matrix_from_metaspace as mfm

metaspace_registry = "2022-08-16_10h14m23s"
api_versions = {'HMDB': ("HMDB", "v4"),
                'LipidMaps': ("LipidMaps", "2017-12-12")}
api_fdr = 0.5
out_dir = "../../../data/whole-ellis/"  # a relative or absolute dir

filename = "whole-ellis"  # no extension

out_file = os.path.join(os.path.expanduser(out_dir),
                        filename)
if not os.path.exists(f'{out_file}.h5ad'):
    print(f"Processing: {metaspace_registry}")
    out_df = mfm.steps_get_united_df(metaspace_registry, api_versions, api_fdr)
    mfm.save_quantifications_anndata(out_df, out_file)
    print(f"Saved: {out_file}.h5ad")
else:
    print("exists, nothing to do")
