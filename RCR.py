import pandas as pd
import networkx as nx
import random
from scipy import special
from collections import defaultdict
import click
import os
from statsmodels.stats import multitest
import csv

"""
Constructs a NetworkX graph.

Input:
    - file : csv file

Output:
    - graph: NetworkX graph
"""
def construct_graph(file):
    df = pd.read_csv(file, sep="\t", header=None, names=["Protein1", "interaction", "Protein2"])
    df_interactions = df.replace("in-complex-with", +1)
    df_interactions = df_interactions.replace("controls-expression-of", -1)
    df_interactions = df_interactions.replace("controls-state-change-of", -1)

    G = nx.DiGraph()

    for i in range(len(df_interactions)):
        prot1 = df_interactions.iloc[i, 0]
    prot2 = df_interactions.iloc[i, 2]
    interaction = df_interactions.iloc[i, 1]
    G.add_node(prot1)
    G.add_node(prot2)
    G.add_edge(prot1, prot2)
    G[prot1][prot2]['relation'] = interaction

    for n, nbrs in G.adj.items():
        for nbr, eattr in nbrs.items():
            relation = eattr['relation']

    return G

"""
Randomly assigns labels of [-1,0,1] to nodes in a graph
Labels:
-1 : Downregulated
0 : No change
+1 : Upregulated

Input:
    - graph : the graph consisting of protein nodes 

Output:
    - prints list of nodes with associated attribute label
"""
def random_node_labels(graph):
    for node in graph.nodes():
        random_label = random.randint(-1,1)
        graph.nodes[node]['label'] = random_label
    print(graph.nodes.data())


"""
Caclulates the shortest path between two nodes.

Input:
    - graph : NetworkX graph
    - source : upstream source node


Output:
    - dictionary of shortest path nodes between source node and all other nodes in graph
"""
def shortest_path(graph, source):
    for target in graph.nodes():
        shortest_paths = nx.shortest_path(graph, source)
    return shortest_paths



"""
Check if node labels of source and target node are the same

Input:
    - graph: NetworkX graph
    - source: source upstream node

Output:
    - list of concordant and non-concordant nodes for the source node
"""


def count_concordance(graph, source):
    same_label = False

    nodes_dic = defaultdict(list)

    for target, path_nodes in shortest_path(graph, source).items():

        # check if node labels of source and target are the same
        if G.nodes[source]['label'] * G.nodes[target]['label'] is 1:
            same_label = True

        # multiply the edge labels
        edge_label = [G[path_nodes[i]][path_nodes[i + 1]]['relation'] * G[path_nodes[i]][path_nodes[i + 1]]['relation']
                      for i in range(len(path_nodes) - 1)]

        # edge_label = 1
        # for i in range(len(path_nodes)-1):
        #    temp_edge_label = G[path_nodes[i]][path_nodes[i+1]]['relation']
        #    edge_label *= temp_edge_label

        # concordant node
        if same_label is True and edge_label is +1:
            graph.nodes[target]['concordance'] = +1
            nodes_dic['concordant'].append(target)

        # non-concordant node
        if same_label is False and edge_label is -1:
            graph.nodes[target]['concordance'] = -1
            nodes_dic['non-concordant'].append(target)

        # no change node
        if G.nodes[source]['label'] is 0 and G.nodes[target]['label'] is 0:
            nodes_dic['no change'].append(target)

    return nodes_dic


"""
Returns a dictionary of the nodes of the graph with their according
     - shortest path nodes
     - concordant nodes
     - non-concordant nodes
     - no change nodes

Input:
    - graph

Output:
    - dictionary of nodes 
"""
def nodes_dictionary(graph):
    dic = {}
    for node in graph.nodes():
        dic[node] = {}

        # concordant, non-concordant and no change nodes
        dic[node] = count_concordance(graph, node)

        # shortest path nodes
        dic[node]['shortest_path'] = list(shortest_path(graph, node).keys()

    return dic


"""
Calculates the concordance for an upstream node with its downstream nodes
Probability of getting at least the number of state changes consistent
with the direction
Input:
    - graph
    - p : probability of achieving a result

Output:
    - dictionary of p-values and corrected p-values for concordance
"""
def calculate_concordance(graph, p):
    concordance_dic = {}

    assert 0 <= p and p <= 1, "p must be within [0,1]"

    if hyp_node not in graph.nodes():
        raise ValueError(f"The node {hyp_node} is not in the graph.")

    for hyp_node in graph.nodes():
        # n is number of trials
        n = len(shortest_path(graph, hyp_node).keys())
        # k is number of successful predictions
        k = len(count_concordance(graph, hyp_node)['concordant'])

        bin_coeff = special.binom(n, k)
        concordance = bin_coeff * (p ** k) * (1 - p) ** (n - k)
        concordance_dic[hyp_node] = {}
        concordance_dic[hyp_node]['concordant'] = concordance

    # correction for multiple testing
    reject, pvals_corrected = multitest.multipletests(concordance_dic.values()),alpha=0.05,method='bonferroni')
    corrected_concordance_dic = {}
    for node, pval in zip(graph.nodes(),pvals_corrected):
        concordance_dic[node]['non-concordant'] = pval

    return concordance_dic


"""
Writes the values for nodes, concordant_nodes, non_concordant_nodes, no_change_nodes, p_val, p_val_corrected to a csv file

Input:
- graph
- csv_output : path for output file 

Output:
- csv file
"""
def write_concordance_csv(graph, csv_output):
    rows = []
    for node in graph.nodes():
        node_dict = {}
        node_dict['node'] = node
        node_dict['concordant_nodes'] = len(count_concordance(graph, node)['concordant'])
        node_dict['non_concordant_nodes'] = len(count_concordance(graph, node)['non-concordant'])
        node_dict['no_change'] = len(count_concordance(graph, node)['no change'])
        node_dict['p_val'] = calculate_concordance(graph)[node]['p_val']
        node_dict['p_val_corrected'] = calculate_concordance(graph)[node]['p_val_corrected']

        rows.append(node_dict)

    df_ = pd.DataFrame(rows)

    df_.to_csv(csv_output)

@click.command()
@click.option('--input', help='The input csv file with nodes and interaction.')
@click.option('--output', help='The output path for the concordance csv file.')


if __name__ == '__main__':

    file = "/Users/sophiakrix/git/MechanismEnrichmentLab/NOTCH1_Intracellular.txt"
    # Create parser
    #parser = argparse.ArgumentParser(description='csv file with Protein1, interaction and Protein2 information')

    # Execute the parse_args() method
    #args = parser.parse_args()
    #file = args.Path

    if not os.path.exists(file):
        raise IOError('The path specified does not exist')

    #file = sys.argv[1]

    # TODO: csv_output argument

    g = construct_graph(file)

    print(f"The concordance of each HYP node of the graph is as follows: ")

    write_concordance_csv(g,csv_output)
