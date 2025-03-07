import os

import zonist.aux.matrix_from_metaspace as mfm


""""
Saves the quantifications to anndata file by fetching their registries from Metaspace. 
The registries are accessed with their date_hour as indicated below. 
The .h5ad files are saved
"""


#metaspace_registry = "2024-06-11_16h33m34s"
api_versions = {'HMDB': ("HMDB", "v4"),
                'KEGG': ("KEGG", "v1")}
api_fdr = 0.5
out_dir = "../../../data/gb-brain/"  # a relative or absolute dir


### NEW
metaspace_registries_and_files = {
    "2025-02-14_15h08m00s": "test635_80-400",
    "2025-02-24_10h42m03s": "test635_200-1K"
}

for metaspace_registry, filename in metaspace_registries_and_files.items():
    out_file = os.path.join(os.path.abspath(out_dir),
                            filename)
    if not os.path.exists(f'{out_file}.h5ad'):
        print(f"Processing: {metaspace_registry}")

        out_df = mfm.steps_get_united_df(metaspace_registry,
                                         api_versions, api_fdr)
        mfm.save_quantifications_anndata(out_df, out_file)
        print(f"Saved: {out_file}.h5ad")
    else:
        print("exists, nothing to do")



