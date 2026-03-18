# expected location of the reference pathways subfolder (within data/)

pathways_files_location_info = {
        "KEGG" : {'file': "KEGG_hsa_pathways_compounds_R110.gmt",
                      'subfolder': "metabolites-mappings",
                      'id_type': 'KEGG'},
        "WIKI" : {'file': "WIKIPATHWAYS_human_multi_id_redundant.mmt",
                          'subfolder': "metabolites-mappings",
                          'id_type': 'moleculeIds'}}


#carefully curated palette (inspired by Okabe–Ito + Paul Tol + Glasbey principles):
custom_palette = [
    '#2F4F4F',  # dark slate gray (strong anchor)
    '#E69F00',  # orange
    '#56B4E9',  # sky blue
    '#009E73',  # green
    '#F0E442',  # yellow
    '#0072B2',  # dark blue
    '#D55E00',  # vermillion
    '#66CC99',  # light green/teal
    '#999999',  # gray
    '#00CED1',  # light cyan
    '#8C564B',  # brown
    '#6A3D9A',  # deep purple
    '#B8860B'  # dark goldenrod
]
