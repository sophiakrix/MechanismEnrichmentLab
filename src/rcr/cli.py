# -*- coding: utf-8 -*-

"""Command line interface."""

import click
import os
import networkx as nx
from src.rcr.rcr_functions import write_concordance_csv, construct_graph
from typing import List
from src.rcr.constants import *
import wrapt

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

output_option = click.option(
    '--output',
    help="Path to output file",
    type=click.Path(file_okay=True, dir_okay=False, exists=False)
)

probability_option = click.option(
    '--probability',
    help="Probability of getting at least the number of state changes consistent with the direction",
    type=float
)

separator_option = click.option(
    '--separator',
    help="Separator for the input file.",
    type=str
)

ppi_columns_option = click.option(
    '--ppicolumns',
    help="Column names of the PPI file",
    type=str
)

dgxp_columns_option = click.option(
    '--dgxpcolumns',
    help="Column names of the DGXP file",
    type=str
)

threshold_option = click.option(
    '--threshold',
    help="The threshold for the fold change of gene expression in the dgxp file to filter by",
    type=float
)


# parent node = click.group(), main
@click.group()
def main():
    """Run rcr."""
    click.echo("Display all RCR commands.")


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
@output_option
@separator_option
@ppi_columns_option
@dgxp_columns_option
@threshold_option
def write_ppi_to_csv(ppi: str=PPI_FILE, dgxp: str=DGXP_FILE, output: str=OUTPUT_FILE, p: float=PROBABILITY, sep: str=SEPARATOR, ppi_columns: str=COLUMNS, dgxp_columns: str=DGXPCOLUMNS, threshold: float=THRESHOLD):
    """ Write the results of concordant nodes and p values to a csv file.
    :param ppi: PPI file
    :param dgxp: DGXP file
    :param output: output path
    :param p: probability
    :param sep: separator
    :param ppi_columns: PPI column names
    :param dgxp_columns: DGXP column names
    :param threshold: threshold for fold-change
    """

    print(f"We want to create a folder here: {PROJECT_DIR}.")
    print(f"The folder is being created here: {OUTPUT_DIR}.")

    # convert input type string into list for DGXP_COLUMNS and PPI_COLUMNS
    if ppi_columns:
        ppi_columns = ppi_columns.split(',')

    if dgxp_columns:
        dgxp_columns = dgxp_columns.split(',')



    # imaging that output is an option alreadz, check if has been given by the user, if so, use it
    # if not
    # Check if user has given this option because it is now not mandatory '(ppi
    #  split and get the last slash out of the path. strip the format .strip('txt)
    #  OUTPUT_from_constnats.format(NAME)

    # if user gives a specific output
    if output:
        output = OUTPUT_FILE.format(output)
    # if user does not give an output use ppi file name for output
    elif output == OUTPUT_FILE:
        name = ppi.split('/')[-1].strip('.txt')
        output = OUTPUT_FILE.format(name)
    # if user does not give an input file use default PPI_FILE and OUTPUT_FILE
    elif ppi == PPI_FILE:
        name = PPI_FILE.split('/')[-1].strip('.txt')
        output = OUTPUT_FILE.format(name)

    graph = construct_graph(ppi, dgxp, sep, ppi_columns, dgxp_columns, threshold)

    write_concordance_csv(graph, output, p)

    click.echo(f"Wrote concordance values to file.")


if __name__ == '__main__':
    main()

