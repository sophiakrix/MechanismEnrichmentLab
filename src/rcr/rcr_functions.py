# -*- coding: utf-8 -*-

"""This module contains all functions used in the RCR command line interface."""

import random
from collections import defaultdict

import networkx as nx
import pandas as pd
from scipy import special
from statsmodels.stats import multitest

from src.rcr.constants import *


def read_ppi_file(ppi_file: str = PPI_FILE, sep: str = SEPARATOR, ppi_columns=COLUMNS) -> pd.DataFrame:
    """Read in a PPI file and replace interaction with activate (+1) or deactivate (-1). 
    
    :param ppi_file: PPI File
    :param sep: separator
    :param ppi_columns: PPI column names 
    :return: dataframe with interactions
    """

    df_ppi = pd.read_csv(ppi_file, sep=sep, header=None)
    # assert that input columns are 3
    df_ppi.columns = ppi_columns

    df_ppi.Protein1 = df_ppi['Protein1'].str.lower()
    df_ppi.Protein2 = df_ppi['Protein2'].str.lower()

    df_interactions = df_ppi.replace("in-complex-with", +1)
    df_interactions = df_interactions.replace("controls-expression-of", -1)
    df_interactions = df_interactions.replace("controls-state-change-of", -1)
    # replace all other interaction expressions with 0
    # df_interactions.loc[(df_interactions.loc[COLUMNS[1]] != +1) & (df_interactions.loc[COLUMNS[1]] != -1)] = 0
    df_interactions[COLUMNS[1]].loc[(df_interactions[COLUMNS[1]] != +1)] = 0
    df_interactions[COLUMNS[1]].loc[(df_interactions[COLUMNS[1]] != -1)] = 0

    return df_interactions


def construct_graph_from_ppi(ppi_file: str = PPI_FILE, sep: str = SEPARATOR, ppi_columns=COLUMNS) -> nx.DiGraph:
    """ Construct a graph from a given PPI file.

    :param ppi_file: PPI file
    :param ppi_columns: column names for PPI file
    :return: NetworkX graph
    """

    df_interactions = read_ppi_file(ppi_file, sep, ppi_columns)

    graph = nx.DiGraph()

    for i in range(len(df_interactions)):
        prot1 = df_interactions.loc[i, COLUMNS[0]]
        prot2 = df_interactions.loc[i, COLUMNS[2]]
        interaction = int(df_interactions.loc[i, COLUMNS[1]])
        # G.add_node(prot1)
        # G.add_node(prot2)

        # edge attributes

        graph.add_edge(prot1, prot2, **{RELATION: interaction})
        # G[prot1][prot2][RELATION] = interaction

    return graph


def filter_dgxp(dgxp_file: str = DGXP_FILE, sep: str = SEPARATOR, dgxp_columns=DGXPCOLUMNS,
                threshold: float = THRESHOLD) -> pd.DataFrame:
    """
    Filter a DGXP file for fold-change significant genes.

    :param dgxp_file: DGXP file
    :param dgxp_columns: column names for DGXP file
    :param threshold: threshold for fold-change
    :return: DataFrame with filtered genes acc. to fold-change
    """
    # node attributes
    df_dgxp = pd.read_csv(dgxp_file, sep=sep, index_col=GENE_GEO_COLUMN)

    #  header=None, names = dgxp_columns ?

    # select only necessary columns
    df_dgxp = df_dgxp[GEO_COLUMNS].copy()
    # change column names accordingly

    # make dictionary from dgxp column names to those, dictionary to constants

    # dictionary for mapping GEO_COLUMNS to DGXPCOLUMNS

    df_dgxp.columns = DGXPCOLUMNS
    df_dgxp = df_dgxp.dropna()
    df_dgxp.index = df_dgxp.index.str.lower()

    # filter dgxp

    # filter nodes with p-value < 0.05
    df_dgxp_pval_filtered = df_dgxp.loc[df_dgxp[DGXPCOLUMNS[1]] < PVALUE]

    # filter nodes with threshold (given by user)
    df_dgxp_thr_filtered = df_dgxp_pval_filtered.loc[abs(df_dgxp_pval_filtered[DGXPCOLUMNS[0]]) > threshold]

    # set fold change labels from float to +1 or -1
    df_dgxp_thr_filtered[FOLD_CHANGE_POSNEG] = df_dgxp_thr_filtered[FOLD_CHANGE].apply(lambda x: +1 if x > 0 else -1)

    #df_dgxp_thr_filtered.loc[df_dgxp_thr_filtered[DGXPCOLUMNS[0]] > 0, FOLD_CHANGE_POSNEG] = +1
    #df_dgxp_thr_filtered.loc[df_dgxp_thr_filtered[DGXPCOLUMNS[0]] < 0, FOLD_CHANGE_POSNEG] = -1

    #df_dgxp_thr_filtered[FOLD_CHANGE].loc[(df_dgxp_thr_filtered[DGXPCOLUMNS[0]] > 0)] = +1
    #df_dgxp_thr_filtered[FOLD_CHANGE].loc[(df_dgxp_thr_filtered[DGXPCOLUMNS[0]] < 0)] = -1

    return df_dgxp_thr_filtered


def create_gene_to_fold_change_dict(dgxp_file: str = DGXP_FILE, ppi_file: str = PPI_FILE, ppi_columns=COLUMNS,
                                    sep: str = SEPARATOR, dgxp_columns=DGXPCOLUMNS,
                                    threshold: float = THRESHOLD) -> dict:
    """Create a dictionary that is mapping the genes in the PPI_FILE to the genes in the DGXP file and take
    key = gene,
    value = the value of the fold-change of the dgxp file
    :param dgxp_file: DGXP file
    :param ppi_file: PPI file
    :param ppi_columns: PPI column names
    :param sep: separator
    :param dgxp_columns: DGXP column names
    :param threshold: threshold for fold change
    :return: dictionary
    """

    df_dgxp = filter_dgxp(dgxp_file, sep, dgxp_columns, threshold)

    df_ppi = read_ppi_file(ppi_file, sep, ppi_columns)

    gene_to_fc_dict = {}

    for node in df_ppi.COLUMNS[0]:
        print('Create gene to fc dict',node)
        if node in df_dgxp.index:
            print(node,' in dg_dgxp index')
            gene_to_fc_dict[node] = {}
            gene_to_fc_dict[node][LABEL] = df_dgxp.loc[node, FOLD_CHANGE_POSNEG]

    return gene_to_fc_dict


def set_node_label(graph: nx.DiGraph, dgxp_file: str = DGXP_FILE, ppi_file: str = PPI_FILE, ppi_columns: str = COLUMNS,
                   sep=SEPARATOR, dgxp_columns=DGXPCOLUMNS,
                   threshold: float = THRESHOLD) -> nx.DiGraph:
    """
    Set the attribute LABEL of nodes according to the fold-change to +1 or -1.

    :param dgxp_file:
    :param sep:
    :param threshold:
    :param dgxp_columns:
    :param graph: NetworkX graph
    :return: graph
    """

    gene_to_fc_dict = create_gene_to_fold_change_dict(dgxp_file, ppi_file, ppi_columns, sep, dgxp_columns, threshold)

    # create a container of tuples to add nodes with label fold change to graph
    all_tup = []

    for node, attr_dict in gene_to_fc_dict.items():
        tup = (node, attr_dict)
        all_tup.append(tup)
        graph.add_nodes_from(nodes_for_adding=all_tup)

    for node in graph.nodes():
        if graph[node][LABEL] is None:
            raise KeyError(f"The node {node} has not been labeled.")

    return graph



def construct_graph(ppi_file: str = PPI_FILE, dgxp_file: str = DGXP_FILE, sep=SEPARATOR, ppi_columns=COLUMNS,
                    dgxp_columns=DGXPCOLUMNS, threshold: float = THRESHOLD) -> nx.DiGraph:
    """
    Construct a NetworkX graph.

    :param ppi_file: tsv file about PPI
    :param dgxp_file: tsv file about differential gene expression
    :param sep:
    :param ppi_columns: list of columns names of PPI file
    :param dgxp_columns: list of column names of DGXP file
    :param threshold:
    :return: NetworkX graph
    """
    # construct the graph
    graph = construct_graph_from_ppi(ppi_file, sep, ppi_columns)

    # label nodes as up- or -downregulated
    node_labeled_graph = set_node_label(graph, dgxp_file, sep, dgxp_columns, threshold)

    return node_labeled_graph


def random_node_labels(graph):
    """
    Randomly assigns labels of [-1,0,1] to nodes in a graph
    Labels:
    -1 : Downregulated
    0 : No change
    +1 : Upregulated

    :param graph: the graph consisting of protein nodes
    """
    for node in graph.nodes():
        random_label = random.randint(-1, 1)
        graph.nodes[node][LABEL] = random_label
    print(graph.nodes.data())


def shortest_path(graph, source) -> dict:
    """
    Caclulates the shortest path between two nodes.

    :param graph: NetworkX graph
    :param source: upstream source node
    :return: dictionary of shortest path nodes between source node and all other nodes in graph
    """
    # for target in graph.nodes():
    shortest_paths = nx.shortest_path(graph, source)
    return shortest_paths


def count_concordance(graph: nx.DiGraph, source) -> dict:
    """
    Check if node labels of source and target node are the same and return dictionary with
    concordant, non-concordant and no change nodes.

    :param graph: NetworkX graph
    :param source: source upstream node
    :return: dictionary of concordant, non-concordant and no change nodes for the source node
    """
    same_label = False

    nodes_dic = defaultdict(list)

    for target, path_nodes in shortest_path(graph, source).items():

        # check if node labels of source and target are the same
        if graph.nodes[source][LABEL] * graph.nodes[target][LABEL] is 1:
            same_label = True

        # multiply the edge labels
        edge_label = [
            graph[path_nodes[i]][path_nodes[i + 1]][RELATION] * graph[path_nodes[i]][path_nodes[i + 1]][RELATION]
            for i in range(len(path_nodes) - 1)]

        # edge_label = 1
        # for i in range(len(path_nodes)-1):
        #    temp_edge_label = G[path_nodes[i]][path_nodes[i+1]]['relation']
        #    edge_label *= temp_edge_label

        # concordant node
        if same_label is True and edge_label is +1:
            graph.nodes[target][CONCORDANCE] = +1
            nodes_dic[CONCORDANT].append(target)

        # non-concordant node
        if same_label is False and edge_label is -1:
            graph.nodes[target][CONCORDANCE] = -1
            nodes_dic[NONCONCORDANT].append(target)

        # no change node
        if graph.nodes[source][LABEL] is 0 and graph.nodes[target][LABEL] is 0:
            nodes_dic[NOCHANGE].append(target)

    return nodes_dic


def nodes_dictionary(graph) -> nx.DiGraph:
    """
    Return a dictionary of the nodes of the graph with their according
         - shortest path nodes
         - concordant nodes
         - non-concordant nodes
         - no change nodes

    :param graph: NetworkX graph
    :return: dictionary of nodes
    """
    dic = {}
    for node in graph.nodes():
        dic[node] = {}

        # concordant, non-concordant and no change nodes
        dic[node] = count_concordance(graph, node)

        # shortest path nodes
        dic[node][SHORTESTPATH] = list(shortest_path(graph, node).keys())

    return dic


def calculate_concordance(graph, p: float = PROBABILITY) -> dict:
    """
    Calculates the concordance for an upstream node with its downstream nodes
    Probability of getting at least the number of state changes consistent
    with the direction

    :param graph: NetworkX graph
    :param p: probability of achieving a result
    :return: dictionary of p-values and corrected p-values for concordance
    """
    concordance_dic = {}

    assert 0 <= p <= 1, "p must be within [0,1]"

    for hyp_node in graph.nodes():
        if hyp_node not in graph.nodes():
            raise ValueError(f"The node {hyp_node} is not in the graph.")
        # n is number of trials - nochange
        all = len(shortest_path(graph, hyp_node).keys())
        no_change = len(count_concordance(graph, hyp_node)[NOCHANGE])
        n = all - no_change
        # k is number of successful predictions
        k = len(count_concordance(graph, hyp_node)[CONCORDANT])

        bin_coeff = special.binom(n, k)
        concordance = bin_coeff * (p ** k) * (1 - p) ** (n - k)
        concordance_dic[hyp_node] = {}
        concordance_dic[hyp_node][PVAL] = concordance

    # correction for multiple testing
    reject, pvals_corrected = multitest.multipletests(concordance_dic.values(), alpha=0.05, method='bonferroni')
    for node, pval in zip(graph.nodes(), pvals_corrected):
        concordance_dic[node][PVALCORRECTED] = pval

    return concordance_dic


def write_concordance_csv(graph, csv_output: str = OUTPUT_FILE, p: float = PROBABILITY):
    """
    Write the values for nodes, concordant_nodes, non_concordant_nodes, no_change_nodes, p_val, p_val_corrected to a csv file

    :param graph: NetworkX graph
    :param csv_output: path for output file
    :param p: probability of achieving a result
    """
    rows = []
    for node in graph.nodes():
        node_dict = {
            NODE: node,
            CONCORDANT: len(count_concordance(graph, node)[CONCORDANT]),
            NONCONCORDANT: len(count_concordance(graph, node)[NONCONCORDANT]),
            NOCHANGE: len(count_concordance(graph, node)[NOCHANGE]),
            PVAL: calculate_concordance(graph)[node][PVAL],
            PVALCORRECTED: calculate_concordance(graph, p)[node][PVALCORRECTED]
        }

        rows.append(node_dict)

    df_ = pd.DataFrame(rows)
    df_.to_csv(csv_output)
