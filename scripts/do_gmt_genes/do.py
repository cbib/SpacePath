import os
import requests
import re
from tqdm import tqdm
import pandas as pd
import sspa as sspa


"""
NOTE:
function adapted from old version v1.0.0 of sspa package ,
because the versions >= 1.2 (and latest) do not support KEGG gene symbols
"""


def local_download_KEGG_genes(organism, filepath=None):
    '''
    credits CWieder sspa v1.0.0 2023, with str split changes by johaGL 2024
    ---------------
    Function for KEGG pathway download FOR GENES
    Args:
        organism (str): KEGG 3 letter organism code
        filepath (str): filepath to save pathway file to, default is None - save to variable

    Returns:
        GMT-like pd.DataFrame containing KEGG pathways
    '''
    print("Beginning KEGG download...")
    # get all pathways
    url = 'http://rest.kegg.jp/list/pathway/' + organism
    # change organism name
    data = requests.get(url)
    pathways = data.text
    pathways = pathways.split("\n")
    pathways = filter(None, pathways)
    pathway_dict = dict()

    for path in pathways:
        path = path.split("\t")
        name = path[1]
        pathid = re.search(r"(.*)", path[0]).group(1)
        pathway_dict[pathid] = name

    # get compounds for each pathway
    base_url = 'http://rest.kegg.jp/get/'

    pathway_ids = [*pathway_dict]
    pathway_names = list(pathway_dict.values())

    # get release details
    release_data = requests.get('http://rest.kegg.jp/info/kegg')
    version_no = release_data.text.split()[9][0:3]

    pathway_features_mapping = dict()

    for index, i in enumerate(tqdm(pathway_ids)):
        myoutputlist = []
        current_url = base_url + "pathway:" + i
        # parse the pathway description page
        page = requests.get(current_url)
        lines = page.text.split("\n")

        try:
            search_start = \
            [lines.index(i) for i in lines if i.startswith("GENE")][0]
            search_field_end = [lines.index(i) for i in lines if
                               i.startswith("COMPOUND") or i.startswith(
                                   "REFERENCE")][0]
            good_lines = lines[search_start:search_field_end]
            prefirst = good_lines.pop(0).split("; ")[0]
            first_feature = prefirst.split(" ")[-1]

            myoutputlist.append(first_feature)
            myoutputlist = myoutputlist + [i.split("; ")[0].split()[-1] for i in good_lines]
            # print(i)
            # print(myoutputlist)
            # print()
            pathway_features_mapping[i] = list(set(myoutputlist))
        except IndexError:
            pathway_features_mapping[i] = []

        # remove empty pathway entries
        pathway_features_mapping = {k: v for k, v in
                                    pathway_features_mapping.items() if v}

    # create GMT style file, genes
    df = pd.DataFrame.from_dict(pathway_features_mapping, orient='index')
    # map pathway names onto first column
    df.insert(0, 'Pathway_name', df.index.map(pathway_dict.get))

    if filepath:
        fpath = filepath + "/KEGG_" + organism + "_pathways_genes_R" + str(
            version_no) + ".gmt"
        df.to_csv(fpath, sep="\t", header=False)
        print("KEGG DB file saved to " + fpath)
    print("Complete!")

    return df


def ancient_fun(organism, filepath):
    pathways_df = local_download_KEGG_genes(organism, filepath)
    print(pathways_df)
    return pathways_df


## start
print("downloading kegg genes per pathway, all pathways")
sspa.process_pathways.process_kegg(
                    organism="hsa", infile=None, download_latest=True,
                    filepath="tmp")
print("downloading kegg metabolites per pathway, all pathways")
df = ancient_fun(organism="hsa", filepath="tmp" )



#  --------------- this FAILED :   -------------------
# from gseapy import Msigdb
# msig = Msigdb()
# gmt = msig.get_gmt(category='mh.all', dbver="2023.2.Mm")
# #  https://docs.gsea-msigdb.org/#MSigDB/Release_Notes/MSigDB_2024.1.Hs/
#
# #u = len(gmt)
# t = 55
# p= 0