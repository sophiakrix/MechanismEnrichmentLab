import unittest

import src.rcr.rcr_functions as rcr
from src.rcr.constants import PPIFILE, DGXPFILE, RELATION, CONCORDANT, NONCONCORDANT, NOCHANGE


class TestRCR(unittest.TestCase):
    """This tests all of the rcr functions."""

    def test_construct_graph_from_ppi(self):
        graph = rcr.construct_graph_from_ppi(PPIFILE)

        nodes = list(graph.nodes())
        edges = list(graph.edges())
        for (u, v, edge_attr) in graph.edges().data(RELATION):
            self.assertIsInstance(
                edge_attr,
                int,
                msg=f'Edge between {u} and {v} is not of type int.')
            self.assertIn(
                edge_attr,
                [-1, 0, +1],
                msg=f'Edge between {u} and {v} is not of in [-1,0,+1].')

        self.assertTrue(
            nodes,
            msg='The graph does not have any nodes.')

        self.assertTrue(
            edges,
            msg='The graph does not have any edges.')

        for node in nodes:
            self.assertIsInstance(
                node,
                str,
                msg=f'The node {node} is not of type str. Please change to type str.')


def test_filter_dgxp(self):
    df = rcr.filter_dgxp(DGXPFILE)

    # check if df not emtpy
    self.assertFalse(df.empty, msg='The dgpx file is empty or it could not be loaded.')

    # check if column 'interaction' has label in [-1,+1]
    for index, row in df.iterrows():
        self.assertIn(
            row[1],
            [-1, +1],
            msg=f'The gene at index {index} does not have a fold-change within [-1, +1].')


def test_set_node_label(self):
    graph = rcr.construct_graph_from_ppi(PPIFILE)

    rcr.set_node_label(graph, DGXPFILE)

    self.assertTrue(
        graph[node][LABEL],
        [-1, +1],
        msg=f'The attribute LABEL of node {node} is not within [-1,+1].')


def test_count_concordance(self):
    graph = rcr.construct_graph_from_ppi(PPIFILE)

    for node in graph.nodes():
        nodes_dic = rcr.count_concordance(graph, node)

        # check if dic is empty
        self.assertTrue(nodes_dic, msg=f'The dictionary with concordance values of node {node} is empty.')

        # check if node label CONCORDANCE is within -1, +1
        self.assertTrue(
            nodes_dic[CONCORDANT],
            +1,
            msg=f'The attribute CONCORDANT of node {node} is not +1.')

        self.assertTrue(
            nodes_dic[NONCONCORDANT],
            -1,
            msg=f'The attribute NONCONCORDANT of node {node} is not -1.')

        self.assertTrue(
            nodes_dic[NOCHANGE],
            0,
            msg=f'The attribute NOCHANGE of node {node} is not 0.')

def test_calculate_concordance(self):
    graph = rcr.construct_graph_from_ppi(PPIFILE)
    



