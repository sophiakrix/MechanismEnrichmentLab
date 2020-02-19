# -*- coding: utf-8 -*-

import os

"""This module contains all the constants used in RCR package."""

RELATION = 'relation'
HERE = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname("__file__"))))
DIRECTORY = os.path.join(HERE,'output')
COLUMNS = ["Protein1", "interaction", "Protein2"]
DGXPCOLUMNS = ["gene", "fold-change", "p-value"]
SEPARATOR = '\t'
LABEL = 'label'
CONCORDANCE = 'concordance'
CONCORDANT = 'concordant'
NONCONCORDANT = 'non-concordant'
NOCHANGE = 'no change'
SHORTESTPATH = 'shortest_path'
NODE = 'node'
PVAL = 'p_val'
PVALCORRECTED = 'p_val_corrected'
PROBABILITY = 0.5
THRESHOLD = 1.2
PPIFILE = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname("__file__")))),'data', 'example.txt')
DGXPFILE = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname("__file__")))),'data', 'example_dgxp.txt')
TESTOUTPUT = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname("__file__"))), 'output', f'{__file__}.csv'))

