# johaGL 2025
# cleaning spots in semimanual fashion for specific tissue (Nature 2025)

import pandas as pd



from spots_preclassif_TD import run_custom_spot_qual

data_centralized_folder = "../../../data"
my_csv_table_path = "../../datasets_parameters/public-datasets-parameters.tsv"
figures_centralized_folder = "../../../figures/"

# start
datasets_parameters_df = pd.read_csv(my_csv_table_path, header=0,
                                     index_col=None, sep='\t')

datasets_parameters_df = datasets_parameters_df.loc[
    datasets_parameters_df['dataset_name'].isin(
        ["Brain-healthy-mouse", "Brain-lmrange-mouse"]), :].copy()  #

weights_spatial_graph = {
    "Brain-healthy-mouse": 0.3,
    "Brain-lmrange-mouse" : 0.3
}
answer = input("continue? y/n")
if answer in ["y", "yes"]:
    run_custom_spot_qual(datasets_parameters_df,
                data_centralized_folder,
                figures_centralized_folder, weights_spatial_graph)

