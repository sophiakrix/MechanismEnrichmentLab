# -*- coding: utf-8 -*-

"""Command line interface."""

import click
import os
import networkx as nx


# TODO : do i need to import rcr_functions.py?



ppi_option = click.option(
    '--ppi',
    help="Path to tab-separated PPI file",
    type=click.Path(file_okay=True, dir_okay=False, exists=True),
    required=True
)

dgxp_option = click.option(
    '--dgxp',
    help="Path to tab-separated differential gene expression file",
    type=click.Path(file_okay=True, dir_okay=False, exists=True),
    required=True
)

output_option = click.option(
    '--output',
    help="Path to csv output file that contains concordance values for HYPs",
    type=click.Path(file_okay=True, dir_okay=False),
    default=HERE
)

# parent node = click.group(), main
@click.group()
def main():
    """Run rcr."""
    print("Display all RCR commands.")


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

    # TODO: do options need to be placed here? --> Vinay's code

# child child node = command of group reverse_causal_reasoning
@reverse_causal_reasoning.command()
@ppi_option
# TODO: add @dgxp_option
def create_graph(ppi: str): -> nx.DiGraph()
    """Create a NetworkX graph from a tsv PPI file."""
    graph = construct_graph


@ppi_option
def write_ppi_to_csv(ppi):
    """ Write the results of concordant nodes and p values to a csv file. """
    HERE = os.path.abspath(os.path.dirname(__file__))
    output = os.path.join(os.getcwd(), ppi+'.csv')

    try:
        os.mkdir(HERE)
    except OSError:
        print("Creation of the directory %s failed" % path)
    else:
        print("Successfully created the directory %s " % path)

    write_concordance_csv(graph, output)


if __name__ == '__main__':
    main()

