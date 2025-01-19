import os

import zonist.aux.matrix_from_metaspace as mfm


""""
Saves the quantifications to anndata file by fetching their registries from Metaspace. 
The registries are accessed with their date_hour as indicated below. 
The .h5ad files are saved with the names 'UBcalibrX' where X is the number of
 the dataset:
 
2024-04-30_16h21m48s : UBcalibr1
2024-08-13_09h45m45s: UBcalibr2
2024-07-18_09h33m42s: UBcalibr3
2024-06-11_16h33m34s: UBcalibr4
2024-04-19_11h59m34s: UBcalibr5
"""


#metaspace_registry = "2024-06-11_16h33m34s"
api_versions = {'HMDB': ("HMDB", "v4"),
                'LipidMaps': ("LipidMaps", "2017-12-12")}
api_fdr = 0.5
out_dir = "../../../data/UBplatform/"  # a relative or absolute dir

#filename = "platform1_mx"  # no extension

### NEW
metaspace_registries_and_files = {
    "2024-04-30_16h21m48s": "UBcalibr1",
    "2024-08-13_09h45m45s": "UBcalibr2",
    "2024-07-18_09h33m42s": "UBcalibr3",
    "2024-06-11_16h33m34s": "UBcalibr4",
    "2024-04-19_11h59m34s": "UBcalibr5"
}

for metaspace_registry, filename in metaspace_registries_and_files.items():
    print(f"Processing: {metaspace_registry}")
    out_file = os.path.join(os.path.expanduser(out_dir),
                            filename)
    out_df = mfm.steps_get_united_df(metaspace_registry,
                                     api_versions, api_fdr)
    mfm.save_quantifications_anndata(out_df, out_file)
    print(f"Saved: {out_file}_mx.h5ad")

