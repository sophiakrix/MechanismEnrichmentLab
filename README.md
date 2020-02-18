# MechanismEnrichmentLab

Input
======================

- file : a csv file with the columns 'Protein1', 'interaction', 'Protein2'

Output
======================

- graph : NetworkX graph

The graph needs to have the following attributes:

Node attributes:
    'label' : [-1,0,+1] the label indicating downregulation, no change and upregulation

Edge attributes:
    'relation' : [-1,+1] indicating a negative or positive regulation of a node to the other
    
