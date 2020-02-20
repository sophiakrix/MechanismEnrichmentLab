# -*- coding: utf-8 -*-

"""Command line interface."""

import click
import os
import networkx as nx
from .rcr_functions import write_concordance_csv, construct_graph
from .constants import RCR_DIR, DIRECTORY
from typing import List
from src.rcr.constants import *

ppi_option = click.option(
    '--ppi',
    help="Path to tab-separated PPI file",
    type=click.Path(file_okay=True, dir_okay=False, exists=True)
)


dgxp_option = click.option(
    '--dgxp',
    help="Path to tab-separated differential gene expression file",
    type=click.Path(file_okay=True, dir_okay=False, exists=True)
)

probability_option = click.option(
    '--probability',
    help="Probability of getting at least the number of state changes consistent with the direction",
    type=click.Path(file_okay=True, dir_okay=False, exists=True),
    required=True
)

separator_option = click.option(
    '--separator',
    help="Separator for the input file.",
    type=str,
)

ppi_columns_option = click.option(
    '--ppicolumns',
    help="Column names of the PPI file",
    type=List[str, str, str]
)

dgxp_columns_option = click.option(
    '--dgxpcolumns',
    help="Column names of the DGXP file",
    type=List[str, str, str]
)

threshold_option = click.option(
    '--threshold',
    help="The threshold for the fold change of gene expression in the dgxp file to filter by",
    type=float,
)

"""output_option = click.option(
    '--output',
    help="Path to csv output file that contains concordance values for HYPs",
    type=click.Path(file_okay=True, dir_okay=False),
    default=HERE
)
"""

# parent node = click.group(), main
@click.group()
def main():
    """Run rcr."""
    click.echo("Display all RCR commands.")


# child node = main.command()
#@ppi_option
#@output_option
#@main.command()
#def runrcr(ppi: str, output: str):
#    """RCR module."""
#    print(ppi)
#    print(output)


# child node = main.group()
@main.group()
def reverse_causal_reasoning():
    """ Does concordance calculation on a graph."""


# child child node = command of group reverse_causal_reasoning
@reverse_causal_reasoning.command()
@ppi_option
# TODO: add @dgxp_option
def create_graph(ppi: str):
    """Create a NetworkX graph from a tsv PPI file."""


@reverse_causal_reasoning.command()
@ppi_option
@dgxp_option
@probability_option
@separator_option
@ppi_columns_option
@dgxp_columns_option
@threshold_option
def write_ppi_to_csv(ppi: str, dgxp: str, p: float, sep: str, ppi_columns: List[str, str, str], dgxp_columns: List[str, str, str], threshold: float):
    """ Write the results of concordant nodes and p values to a csv file. """

    # TODO take from constants
    print(f"We want to create a folder here: {PRO}.")
    print(f"The folder is being created here: {directory}.")



    # imaging that output is an option alreadz, check if has been given by the user, if so, use it
    # if not
    # Check if user has given this option because it is now not mandatory '(ppi
    #  split and get the last slash out of the path. strip the format .strip('txt)
    #  OUTPUT_from_constnats.format(NAME)

    graph = construct_graph(ppi, dgxp, sep, ppi_columns, dgxp_columns, threshold)

    write_concordance_csv(graph, output, p)

    click.echo(f"Wrote concordance values to file.")


if __name__ == '__main__':
    main()

