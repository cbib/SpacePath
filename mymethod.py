import sys
import os
import logging
import click
import warnings


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('data_dir', type=click.Path(file_okay=False, dir_okay=True))
@click.argument('figures_dir',
                type=click.Path(file_okay=False, dir_okay=True))


def main(input_file, data_dir, figures_dir):
    """
    \b
    'mymethod': Spatial Metabolomics Exploration Based on Pathway and Families Activity Profiles

    Positional arguments:

    \b
    INPUT_FILE: Relative or absolute path to the table in tab delimited format,
     where each dataset and user-defined parameters for its analysis are stored.

    DATA_DIR: Relative or absolute path, to the folder where the datasets are stored, e.g. data/.

    FIGURES_DIR: Relative or absolute path to output figures directory.
    """
    click.echo(f'Input file: {input_file}')
    click.echo(f'Data directory: {data_dir}')
    click.echo(f'Figures directory: {figures_dir}')

    # start
    try:
        if (input_file is not None) and (
            data_dir is not None
        ) and (figures_dir is not None):
            warnings.simplefilter('ignore', category=FutureWarning)
            from src.tissues_dispatcher import process_tissues
            process_tissues(input_file, data_dir, figures_dir)
    except Exception as e:
        logging.error(e)


if __name__ == '__main__':
    main()
