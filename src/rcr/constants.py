# -*- coding: utf-8 -*-

"""This module contains all the constants used in RCR package."""

import os

RELATION = 'relation'
RCR_DIR = os.path.abspath(os.path.dirname(__file__))
SRC_DIR = os.path.abspath(os.path.join(RCR_DIR, os.path.pardir))
PROJECT_DIR = os.path.abspath(os.path.join(SRC_DIR, os.path.pardir))

DATA_DIR = os.path.join(PROJECT_DIR, 'data')

PPI_FILE = os.path.join(DATA_DIR, 'example.txt')
DGXP_FILE = os.path.join(DATA_DIR, 'example_dgxp.txt')

OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")

os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_FILE = os.path.join(OUTPUT_DIR, '{}_results.csv')

INDEX = 'index'
COLUMNS = ['Protein1', "interaction", 'Protein2']
DGXPCOLUMNS = ['fold-change', 'p-value']
GENE_GEO_COLUMN = 'Gene.symbol'
GEO_COLUMNS = ['logFC', 'adj.P.Val']
FOLD_CHANGE_POSNEG = 'fold-change pos-neg'
GENE = 'gene'
FOLD_CHANGE = 'fold-change'
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
THRESHOLD = 0.2
PVALUE = 1.0

