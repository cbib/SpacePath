
from typing import List, Dict

#
mz_decimal_places = 6


supported_annotations: List[str] = ["LipidMaps",
                                    "LIPID_MAPS", # ancient db version too
                                    "CoreMetabolome",  # new
                                    "HMDB",
                                    "KEGG",
                                    "custom_kegg"]

columns_ref_map_expected: List[str] = ["initial_id_unconverted", "KEGG", "database"]


not_metabolic_reactions_pathways = [
    "autophagy - other",  "amoebiasis", "cgmp-pkg signaling",
    "phagocytosis",  "fc gamma", "ri signaling", "mapk",
    "wnt", "notch", "hedgehog",  "erbb", "tgf", "hippo", "vegf", "apelin",
    "jak-stat", "kappa", "tnf", "rap1",  "foxo",
    "rage signal", "oxytocin", "relaxin", "prolactin", "renin","cytokine",
     "gnrh", "salivary secretion", "acid secretion",
    "insulin secretion", "parathyroid", "aldosterone" ,"leukemia",
    "syndrome", "thyroid", "endocrine",  "cushing", "ovarian", "contraction",
    "atherosclerosis", "cardiomyo", "arthritis", "inflammatory", "platelet",
    "lupus erythematosus",  "kaposi", "asthma", "disease", "depression",
    "longevity regulating", "circadian", "thermogenesis", "phototransduction",
    "trypanosomiasis", "shigellosis", "infection", "diabetes",
    "leishmania", "toxoplasmosis", "legionel", "malaria", "tuberculosis",
    "coccus", "virus", "virion", "viral", "measles", "dna adducts",
    "meiosis", "checkpoint", "resistance", "cell",
    "pathways in cancer", "morphine addiction", "gentamicin",  # new
    "xenobiotics", "mineral absorption", "digestion and absorption", # new
    "glucagon signaling", "taste transduction", "drug" # new
]
# Explanation of 'not_metabolic_reactions_pathways': any pathway name
# that contains the specified substrings above, is classified as being
# NOT strictly associated to the sets of metabolic reactions.
# Such substrings refer to disease, or generalistic sets or terms.


supported_mappings: List[str] = ["LipidMaps",
                                 "CoreMetabolome",  # new
                                 "HMDB"]  # this list >
# > states the types of ref tables that match all not-kegg IDs with KEGG ones