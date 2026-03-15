from typing import List


mz_decimal_places = 6

supported_annotations: List[str] = ["LipidMaps",
                                    "LIPID_MAPS",  # ancient db version too
                                    "CoreMetabolome",
                                    "HMDB",
                                    "Swisslipids",
                                    "KEGG",
                                    "custom_kegg"]

columns_kegg_conversion_expected: List[str] = [
    "initial_id_unconverted", "KEGG", "database"]

not_metabolic_reactions_pathways = [
    "autophagy - other", "amoebiasis", "cgmp-pkg signaling",
    "phagocytosis", "fc gamma", "ri signaling", "mapk",
    "wnt", "notch", "hedgehog", "erbb", "tgf", "hippo", "apelin",
    "jak-stat", "kappa", "tnf", "rap1", "foxo",
    "rage signal", "oxytocin", "relaxin", "prolactin", "renin", "cytokine",
    "gnrh", "salivary secretion", "acid secretion",
    "insulin secretion", "parathyroid", "aldosterone", "leukemia",
    "syndrome", "thyroid", "endocrine", "cushing", "ovarian", "contraction",
    "atherosclerosis", "cardiomyo", "arthritis", "inflammatory", "platelet",
    "lupus erythematosus", "kaposi", "asthma", "disease", "depression",
    "longevity regulating", "circadian", "thermogenesis", "phototransduction",
    "trypanosomiasis", "shigellosis", "infection", "diabetes",
    "leishmania", "toxoplasmosis", "legionel", "malaria", "tuberculosis",
    "coccus", "virus", "virion", "viral", "measles", "dna adducts",
    "meiosis", "checkpoint", "resistance", " cell ",
    "hepatocell", "cells", "pathways in cancer",
    "addiction", "gentamicin", "xenobiotics",
    "mineral absorption", "digestion and absorption",
    "glucagon signaling", "taste transduction", "drug",
    "lectin", "recombination", "vitamins and cofactors", "disorder",
    "number variation", "hypothe", "senescence", "mitotic", "mitosis",
    "rhodopsin", "nf-kb", "spina bifida", "nucleotide gpcr", "telomere",
    "follicle devel", "functionalization of compounds", "visual cycle",
    "alcoholism", "hepatitis", "pertussis", "sars-cov-2",
    "autism", "signaling events", "aspirin", "cyclosporin", "sulindac",
    "rotenone", "overview", "workshop"
]  # spaces comprised for each string
# Explanation of 'not_metabolic_reactions_pathways': any pathway name
# that contains the specified substrings above, is classified as being
# NOT strictly associated to the sets of metabolic reactions, or too vague.
# Such substrings refer to disease, or generalistic sets or terms.


supported_mappings: List[str] = ["LipidMaps",
                                 "CoreMetabolome",
                                 "HMDB"]  # this list >
# > states the types of ref tables that match all not-kegg IDs with KEGG ones

message_multimz = (
    "Identified cases where one compound ID maps to multiple m/z: "
    "Please run '<metabo_annotations>.keep_best_mz_id_pairs(anndata_mz)'.")
