# Entry point for the CLI
import logging
import warnings

import click

warnings.simplefilter('ignore', category=FutureWarning)


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('data_dir', type=click.Path(file_okay=False,
                                            dir_okay=True))
@click.argument('figures_dir',
                type=click.Path(file_okay=False, dir_okay=True))
@click.option('--mode', default="all")
def main_cli(input_file: str, data_dir: str, figures_dir: str, mode: str):
    """
    \b
    'SpacePath': Spatial Metabolomics Exploration Based on Pathways
     and Families Profiles

    \b
    Positional arguments:

      INPUT_FILE: Relative or absolute path to the table in tab delimited
      format, where each dataset and user-defined parameters for its
      analysis are stored.

      DATA_DIR: Relative or absolute path, to the folder where the datasets
      are stored, e.g. data/.

      FIGURES_DIR: Relative or absolute path to output figures directory.

    \b
    Optional arguments:
      --mode : Optional, way the method is to be run (all steps, or the
      selected step)
    """
    click.echo(f'Input file: {input_file}')
    click.echo(f'Data directory: {data_dir}')
    click.echo(f'Figures directory: {figures_dir}')
    click.echo(f'  --mode: {mode}')

    try:
        if (input_file is None) or (data_dir is None) or (
                figures_dir is None):
            raise ValueError(
                "Missing required arguments: input_file, data_dir, "
                "or figures_dir.")

        from SpacePath.src.orchestrator import DatasetsProcessor

        multiple_analyses = DatasetsProcessor(
            input_file=input_file, data_dir=data_dir,
            figures_dir=figures_dir, mode=mode
        )
        multiple_analyses.process_tissues()

    except FileNotFoundError as e:
        logging.error(f"File not found: {e}")
        raise  # Re-raise to exit with error code
    except NotADirectoryError as e:
        logging.error(f"Directory not found: {e}")
        raise
    except ValueError as e:
        logging.error(f"Invalid input: {e}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        raise


if __name__ == "__main__":
    main_cli()
