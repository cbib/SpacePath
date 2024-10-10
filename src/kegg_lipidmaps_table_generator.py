# johaGL 2024
# usage:
#   cd $MYWDIR
#   outdir=data/db_tables/  # must exist in $MYWDIR
#   python3 $TOOLDIR/src/kegg_lipidmaps_table_generator.py $outdir
import os
import requests
import time
import sys

"""
Obtain a table with columns KEGG, LIPIDMAPS, using KEGG API.
Instead running this script, we advise using the table that we provide, 
a prompt allows the user to interrupt the operation.
Otherwise, will require ~7 hours runtime (queries all KEGG compounds)
"""

def kegg_text_to_compound_list(all_compounds_data: str):
    all_compounds_l = all_compounds_data.replace('\t', '  ').split('\n')

    return [li.split("  ")[0].strip() for li in all_compounds_l]


def retrieve_lipidmaps_identifier(one_metabolite_plain_text: str):
    """ENTRY       C00550                      Compound
    NAME        Sphingomyelin
    FORMULA     C24H49N2O6PR
    REMARK      ATC code: A05BA10
    COMMENT     Generic compound in reaction hierarchy
    REACTION    R01891 R02541 R02542 R02543 R08969
    PATHWAY     map00600  Sphingolipid metabolism

    NETWORK     nt06014  Sphingolipid degradation
    ENZYME      2.3.1.24        2.7.8.3         2.7.8.27        3.1.4.12
                3.1.4.41
    BRITE       Compounds with biological roles [BR:br08001]
                ...
    DBLINKS     PubChem: 3831
                ChEBI: 17636 184744 185051 190464 62490
                LIPIDMAPS: LMSP03010000
    ATOM        34
    """
    txt_whole = one_metabolite_plain_text.split('\n')
    lipidmaps_id_text = ''
    activate_lines_below = False
    lipidmaps_not_found = True
    further_continue = True
    k = 0
    while lipidmaps_not_found and further_continue:
        if txt_whole[k].startswith("DBLINKS"):
            tmp = txt_whole[k].split("    ")[-1] # e.g. :   DBLINKS     PubChem: 3831
            if tmp.startswith("LIPIDMAPS"):
                lipidmaps_id_text = tmp.split(":")[-1].strip()
                lipidmaps_not_found = False
            else:
                activate_lines_below = True
        if activate_lines_below:
            # e.g.:    LIPIDMAPS: LMSP03010000
            tmp = txt_whole[k].split("    ")[-1]
            if tmp.startswith("LIPIDMAPS"):
                lipidmaps_id_text = tmp.split(":")[-1].strip()
                lipidmaps_not_found = False
        if txt_whole[k].startswith("ATOM") or txt_whole[k].startswith("///"):
            further_continue = False
        k += 1

    return lipidmaps_id_text


def write_table_to_file(compounds_ids_l, my_output_file: str):
    txt_out = "KEGG\tLIPIDMAPS\thowmany\n"
    f = open(my_output_file, "a")
    f.write(txt_out)
    for curr_compound in compounds_ids_l:
        one_compound_url = "http://rest.kegg.jp/get/compound:" + curr_compound
        try:
            one_compound_response = requests.get(one_compound_url)
            lipidmapsid_text = retrieve_lipidmaps_identifier(
                one_compound_response.text)
        except Exception as e:
            print(e)
            lipidmapsid_text = ''
        if lipidmapsid_text == '':
            lipidmapsid_nb = 0
        else:
            lipidmapsid_nb = len(lipidmapsid_text.split(' '))
        f.write(f"{curr_compound}\t{lipidmapsid_text}\t{lipidmapsid_nb}\n")
        print(curr_compound)

    f.close()
    print("end")

def run_lipidmaps_table_generator(output_dir):
    print(output_dir)
    print(os.listdir(output_dir))
    print("This script will require ~7 hours (queries all KEGG compounds), "
          "we advise using the table that we provide")
    choice = input("Continue anyway?[y/n]: ")
    if choice.lower() == 'y':
        print("Running kegg to lipidmaps table generator")
        print("Retrieve all compounds")
        base_url = 'http://rest.kegg.jp/list/compound'
        all_compounds_response = requests.get(base_url)
        compounds_ids_l = kegg_text_to_compound_list(
            all_compounds_response.text)

        print("Retrieve lipidmaps id for each of the ",
              len(compounds_ids_l), "compounds")

        out_path = os.path_join(output_dir,
                                f"table_kegg2lipidmaps_{time.time()}.txt")
        write_table_to_file(compounds_ids_l, out_path)
        print("done, saved table to:", out_path)
    elif choice.lower() == 'n':
        raise SystemExit("aborted, exiting")
    elif choice.lower() not in ['y', 'n']:
        raise SystemExit("aborted, exiting")


if __name__ == "__main__" :
    run_lipidmaps_table_generator(sys.argv[1])

