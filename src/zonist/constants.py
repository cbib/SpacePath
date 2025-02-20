
from typing import List, Dict

#

supported_mappings: List[str] = ["LipidMaps", "HMDB"]  # this list >
# > states the types of ref tables that match all not-kegg IDs with KEGG ones

supported_annotations: List[str] = ["LipidMaps", "HMDB", "KEGG", "custom_kegg"]

columns_ref_map_expected: List[str] = ["initial_id_unconverted", "KEGG", "database"]


not_metabolic_reactions_pathways = [
    "autophagy", "efferocytosis", "amoebiasis", "cgmp-pkg signaling",
    "phagocytosis",  "necroptosis", "fc gamma", "ri signaling", "mapk",
    "wnt", "notch", "hedgehog",  "erbb", "tgf", "hippo", "vegf", "apelin",
    "jak-stat", "kappa", "tnf", "ras s", "rap1",  "foxo", "mtor", "pi3k-akt",
    "rage signal", "oxytocin", "relaxin", "estrog", "prolactin", "renin",
    "cytokine", "d signaling", "vegf", "gnrh", "secretion", "leukemia",
    "syndrome", "thyroid", "endocrine",  "cushing", "ovarian", "contraction",
    "atherosclerosis", "cardiomyo", "arthritis", "inflammatory", "platelet",
    "lupus erythematosus",  "kaposi", "asthma", "disease", "depression",
    "longevity regulating", "circadian", "thermogenesis", "phototransduction",
    "trypanosomiasis", "shigellosis", "infection", "immune", "diabetes",
    "leishmania", "toxoplasmosis", "legionel", "malaria", "tuberculosis",
    "coccus", "virus", "virion", "viral", "measles", "dna adducts",
    "meiosis", "checkpoint", "receptor", "resistance", "cell",
    "pathways in cancer"  # new
]
# Explanation of 'not_metabolic_reactions_pathways': any pathway name
# that contains the specified substrings above, is classified as being
# NOT strictly associated to the sets of metabolic reactions.
# Such substrings refer to disease, signaling, regulation or endocrine system.