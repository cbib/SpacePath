# this is a simple auxiliary script for changing the color of clusters plots
import os
import click
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from SpacePath.src.constants import options_palette_supported

from SpacePath.src.helpers import plot_clustering_results, safe_load_processed_data, safe_insert_clusters_column

mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Arial"
plt.rcParams.update({'font.size': 20})


def run_recolorizer(input_file, data_dir, figures_dir, use_palette):
    from SpacePath.src.orchestrator import DatasetsProcessor as MultiSP

    data_centralized_folder = data_dir
    figures_centralized_folder = figures_dir

    datasets_parameters_df = pd.read_csv(input_file,
                                         header=0, index_col=None, sep='\t')


    for i, row in datasets_parameters_df.iterrows():
        # each row represents a dataset and its parameters
        (SAMPLE, annotations_dict, offSample_cutoff, msm_cutoff, fdr_cutoff,
         moran_I_cutoff, moran_padj_cutoff, clustering_parameters,
         pathways_types, min_pathway_size, cutoff_tree,
         min_family_size
         ) = MultiSP.extract_params_per_dataset(row)
        subfolder = SAMPLE
        adata_mz = safe_load_processed_data(SAMPLE, subfolder,
                                             data_centralized_folder )

        str_out = f"{clustering_parameters['choice']}-{clustering_parameters['numeric_parameter']}"
        clusters_df = pd.read_csv(
            os.path.join(data_centralized_folder,
                         subfolder,f'{SAMPLE}-clusters-{str_out}.txt'), header=0)
        adata_mz = safe_insert_clusters_column(adata_mz, clusters_df,
                                               column_interest=str_out)
        plot_folder = os.path.join(figures_centralized_folder, subfolder)

        plot_clustering_results(
            adata_mz, clustering_parameters, SAMPLE,
            plot_folder,
            colors_spots=options_palette_supported[use_palette])


###############################
# cli interface for this script
###############################

@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('data_dir', type=click.Path(file_okay=False, dir_okay=True))
@click.argument('figures_dir',
                type=click.Path(file_okay=False, dir_okay=True))
@click.option('--use_palette', default="default")

def main(input_file, data_dir, figures_dir, use_palette):
    """
    \b
    Re-colorize the clusters (overwrites existing files):

    \b
    Positional arguments:

      INPUT_FILE: Relative or absolute path to the table in tab delimited format,
        where each dataset and user-defined parameters for its analysis are stored.

      DATA_DIR: Relative or absolute path, to the folder where the datasets are stored, e.g. data/.

      FIGURES_DIR: Relative or absolute path to output figures directory.

      use_palette: one of the options:  "default" or "colorblindness"


    """
    click.echo(f'Input file: {input_file}')
    click.echo(f'Data directory: {data_dir}')
    click.echo(f'Figures directory: {figures_dir}')
    click.echo(f'  --use_palette: {use_palette}')
    # start
    try:
        if (input_file is not None) and (
                data_dir is not None
        ) and (figures_dir is not None):
            run_recolorizer(
                input_file, data_dir, figures_dir, use_palette=use_palette)

    except Exception as e:
            print(e)


if __name__ == '__main__':
    main()