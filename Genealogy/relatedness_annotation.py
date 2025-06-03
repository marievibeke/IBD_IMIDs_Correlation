#!/usr/bin/python

import pandas as pd
import networkx as nx 
import gzip

#Load graph in nx - both directed and un-directed
edge_list = pd.read_csv('/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/edge_file_IMID.csv', header = 0, sep=' ')
G_directed = nx.from_pandas_edgelist(edge_list, source = 'ID', target = 'parentID',create_using=nx.DiGraph())
#Note: the direction is wrong, but enables bfs_predecessors

#Find number of generations:
#n_generations = nx.topological_generations(G_directed)
#n_gen = 0
#for generation in n_generations:
#  n_gen +=1
#print('Number of generations: '+ str(n_gen))

#Load in file with pairs 
res = pd.read_table('/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/relatedness_coeffs_all_no_duplicates_IMID.csv.gz', sep = ' ')

#Make file for results
f = gzip.open('/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/relatedness_annotation_IMID.txt.gz', 'wt')
f.write('id1 id2 n_nodes_in_path depth_1 depth_2 R n_paths parent_rel annot i\n') #Header for file -> the i is added to keep track of progress

#Neccesary functions
def table(seq):
  seq = [x for x in seq]
  keys = set(seq)
  return dict([(k, seq.count(k)) for k in keys])

def path_downward(start, end): #Function to find a list of paths from ancestor (start) to descendant (end)
  return [p for p in nx.all_shortest_paths(G_directed, source = end, target = start)]
  
def paths_through_ancestor(ind1, ind2, ancestor):
  identified_paths = []
  ups = path_downward(ancestor, ind1)
  downs = path_downward(ancestor, ind2)
  for u in ups:
    for d in downs:
      if len(u)>1 and len(d)>1 and u[-2] == d[-2]:
        continue
      path = u[:-1]+d[::-1]
      pt = table(path)
      if [pt[k] for k in pt if pt[k]>1]:
        continue
      identified_paths.append(path)
  return identified_paths

def paths(ind1, ind2, common_anc):
  identified_paths = []
  for ancestor in common_anc:
    identified_paths.extend(paths_through_ancestor(ind1, ind2, ancestor))
  return identified_paths

    

#Loop through rows in the table with pairs of individuals
i=0
for index, row in res.iterrows():
  i+=1
  #Find closest common ancestor of the pair
  preds_1 = nx.bfs_predecessors(G_directed, row['id1'])
  preds_1 = [row['id1']]+ [n[0] for n in preds_1]
  preds_2 = nx.bfs_predecessors(G_directed, row['id2'])
  preds_2 = [row['id2']]+ [n[0] for n in preds_2]
  common_preds = set(preds_1).intersection(set(preds_2))
  lca = min(list(common_preds), key=lambda n: preds_1.index(n)) #Lowest common ancestor
  hca = max(list(common_preds), key=lambda n: preds_1.index(n)) #Highest common ancestor
    
  depth_1 = len(nx.shortest_path(G_directed, source =row['id1'], target = lca))
  depth_2 = len(nx.shortest_path(G_directed, source =row['id2'], target = lca)) 
  
  n_nodes_in_path = depth_1+depth_2-1
  r = row['r']
  
  n_paths = len(paths(row['id1'], row['id2'], list(common_preds)))
    
  #Detect number of parent pairs with common ancestors:
  par1 = nx.bfs_predecessors(G_directed, source = row['id1'], depth_limit = 1)
  par1 = [n[0] for n in par1]
  par2 = nx.bfs_predecessors(G_directed, source = row['id2'], depth_limit = 1)
  par2 = [n[0] for n in par2]
  if len(par1) == 0 or len(par2) == 0:
    parent_anc = 0
  else: #You either have 0 or 2 parent links
    par1_1 = par1[0]
    par1_2 = par1[1]
    par2_1 = par2[0]
    par2_2 = par2[1]
    preds1_1 = nx.bfs_predecessors(G_directed, par1_1)
    preds1_1 = [par1_1]+ [n[0] for n in preds1_1]
    preds1_2 = nx.bfs_predecessors(G_directed, par1_2)
    preds1_2 = [par1_2]+ [n[0] for n in preds1_2]
    preds2_1 = nx.bfs_predecessors(G_directed, par2_1)
    preds2_1 = [par2_1]+ [n[0] for n in preds2_1]
    preds2_2 = nx.bfs_predecessors(G_directed, par2_2)
    preds2_2 = [par2_2]+ [n[0] for n in preds2_2]
    anc_list = [len(set(preds1_1).intersection(set(preds2_1))), len(set(preds1_1).intersection(set(preds2_2))), len(set(preds1_2).intersection(set(preds2_1))), len(set(preds1_2).intersection(set(preds2_2)))]
    parent_anc = len([1 for i in anc_list if i>0])
      
  #Annotation of relationship
  annot = 'NA'
  if n_nodes_in_path == 2:
    if round(r,1) == 0.5 and max(depth_1,depth_2) ==2 and min(depth_1, depth_2)==1 and (parent_anc == 2 or parent_anc == 0) and n_paths == 1:
      annot = 'PO'
  elif n_nodes_in_path == 3:
    if round(r, 1) == 0.5 and depth_1 == 2 and depth_2 == 2 and parent_anc == 2 and n_paths == 2:
      annot = 'FS'
    elif round(r, 2) == 0.25 and n_paths == 1:
      if depth_1 == 2 and depth_2 == 2 and parent_anc == 1:
        annot = 'HS'
      elif max(depth_1,depth_2) ==3 and min(depth_1, depth_2)==1 and (parent_anc == 2 or parent_anc == 0):
        annot = '1G'
  elif n_nodes_in_path == 4:
    if round(r, 2) == 0.25 and max(depth_1,depth_2) ==3 and min(depth_1, depth_2)==2 and parent_anc == 2 and n_paths == 2:
      annot = 'Av'
    elif round(r, 3) == 0.125 and n_paths == 1:
      if max(depth_1,depth_2) ==3 and min(depth_1, depth_2)==2 and parent_anc == 1:
        annot = 'HAv'
      elif max(depth_1,depth_2) ==4 and min(depth_1, depth_2)==1 and (parent_anc == 2 or parent_anc == 0):
        annot = '2G'
  elif n_nodes_in_path == 5:
    if round(r, 3) == 0.125 and n_paths == 2:
      if max(depth_1,depth_2) ==4 and min(depth_1, depth_2)==2 and parent_anc == 2:
        annot = '1GAv'
      elif depth_1 == 3 and depth_2 == 3: 
        if parent_anc == 1:
          annot = '1C'
        elif parent_anc == 2:
          annot = 'H1Cx2'
    elif round(r, 4) == 0.0625 and n_paths == 1:
      if max(depth_1,depth_2) ==4 and min(depth_1, depth_2)==2 and parent_anc == 1:
        annot = '1GHAv'
      elif max(depth_1,depth_2) ==5 and min(depth_1, depth_2)==1 and (parent_anc == 2 or parent_anc == 0):
        annot = '3G'
      elif depth_1 == 3 and depth_2 == 3 and parent_anc == 1:
        annot = 'H1C'
  elif n_nodes_in_path == 6:
    if round(r, 4) == 0.0625 and n_paths == 2:
      if max(depth_1,depth_2) ==4 and min(depth_1, depth_2)==3:
        if parent_anc == 1:
          annot = '1C1R'
        elif parent_anc == 2:
          annot = 'H1C1Rx2'
      elif max(depth_1,depth_2) ==5 and min(depth_1, depth_2)==2 and parent_anc == 2:
        annot = '2GAv'
    elif round(r, 5) == 0.03125 and n_paths == 1:
      if max(depth_1,depth_2) ==6 and min(depth_1, depth_2)==1 and (parent_anc == 2 or parent_anc == 0):
        annot = '4G'
      elif max(depth_1,depth_2) ==4 and min(depth_1, depth_2)==3 and parent_anc == 1:
        annot = 'H1C1R'
      elif max(depth_1,depth_2) ==5 and min(depth_1, depth_2)==2 and parent_anc == 1:
        annot = '2GHAv'
  elif n_nodes_in_path == 7:
    if round(r, 5) == 0.03125 and n_paths == 2:
      if depth_1 == 4 and depth_2 == 4:
        if parent_anc == 1:
          annot = '2C'
        elif parent_anc == 2:
          annot = 'H2Cx2'
      elif max(depth_1,depth_2) ==5 and min(depth_1, depth_2)==3:
        if parent_anc == 1:
          annot = '1C2R'
        elif parent_anc == 2:
          annot = 'H1C2Rx2'
      elif max(depth_1,depth_2) ==6 and min(depth_1, depth_2)==2 and parent_anc == 2:
        annot = '3GAv'
    elif round(r, 6) == 0.015625 and n_paths == 1:
      if depth_1 == 4 and depth_2 == 4 and parent_anc == 1:
        annot = 'H2C'
      elif max(depth_1,depth_2) ==5 and min(depth_1, depth_2)==3 and parent_anc == 1:
        annot = 'H1C2R'
  elif n_nodes_in_path == 8:
    if round(r, 6) == 0.015625 and n_paths == 2:
      if max(depth_1,depth_2) ==5 and min(depth_1, depth_2)==4:
        if parent_anc == 1:
          annot = '2C1R'
        elif parent_anc == 2:
          annot = 'H2C1Rx2'
      elif max(depth_1,depth_2) ==6 and min(depth_1, depth_2)==3:
        if parent_anc == 1:
          annot = '1C3R'
        elif parent_anc == 2:
          annot = 'H1C3Rx2'
    elif round(r, 7) == 0.0078125 and max(depth_1,depth_2) ==5 and min(depth_1, depth_2)==4 and parent_anc == 1 and n_paths == 1:
      annot = 'H2C1R'
  elif n_nodes_in_path == 9:
    if round(r, 7) == 0.0078125 and n_paths == 2:
      if max(depth_1,depth_2) ==6 and min(depth_1, depth_2)==4:
        if parent_anc == 1:
          annot = '2C2R'
        elif parent_anc == 2:
          annot = 'H2C2Rx2'
      elif depth_1 == 5 and depth_2 == 5:
        if parent_anc == 1:
          annot = '3C'
        elif parent_anc == 2:
          annot = 'H3Cx2'
  elif n_nodes_in_path == 10 and round(r, 8) == 0.00390625 and max(depth_1,depth_2) ==6 and min(depth_1, depth_2)==5 and n_paths == 2:
    if parent_anc == 1:
      annot = '3C1R'
    elif parent_anc ==2:
      annot = 'H3C1Rx2'
  elif n_nodes_in_path == 11 and round(r, 9) == 0.001953125 and depth_1 == 6 and depth_2 ==6 and n_paths == 2:
    if parent_anc == 1:
      annot = '4C'
    elif parent_anc == 2:
      annot = 'H4Cx2'
  
  #Write to file
  f.write(str(row['id1'])+' '+str(row['id2'])+' '+str(n_nodes_in_path)+' '+str(depth_1)+' '+str(depth_2)+' '+str(r)+' '+str(n_paths)+' '+str(parent_anc)+' '+str(annot)+' '+str(i)+'\n') 
  if i%100000 == 0:
    print('Status: '+str(i))
    f.flush()
      

f.close()
