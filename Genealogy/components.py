#!/usr/bin/python

import networkx as nx 
import pandas as pd

#read the edge list
G = nx.read_edgelist('/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/edge_file_IMID.csv')

#find the components of the graph
components = nx.connected_components(G)

#register id's in order of appearance in the components
indid = []

#register component size
sizes = [] 

for component in components:
    indid.extend(component)
    sizes.append(len(component))

#register famid's according to component size
famid = [] 
for i in range(len(sizes)):
    famid.extend([i+1]*sizes[i]) 

ids = pd.DataFrame({'famid':famid, 'indid':indid})
export_csv = ids.to_csv('/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/ids_IMID.csv', sep=' ', index=None, header=True)
