#!/usr/bin/python

import pandas as pd
import PyAGH
from itertools import compress
import gzip
import networkx as nx

ped = pd.read_csv('/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/ids_file_preped_IMID.csv', header = 0)

f = gzip.open('/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/relatedness_coeffs_less_IMID.txt.qz', 'wt')

for family in ped.iloc[:,0].unique():
  if family != 2:
    print(family)
    #Subset dataset by families, and remove the column:
    ped_sub = ped.loc[ped['famID'] == family]
    ped_sub = ped_sub.drop('famID', axis = 1)
    ped_sub = pd.DataFrame(data=ped_sub)
    
    #Sort the ped
    ped_sub_sort = PyAGH.sortPed(ped_sub)
    
    #Select smaller pedigrees based on founders and their descendants
    edge_list = pd.concat([ped_sub.iloc[:, [0,1]].astype(str), ped_sub.iloc[:, [0,2]].astype(str).rename(columns={"DAM": "SIRE"})])
    edge_list = edge_list.iloc[:, [0,1]]

    G = nx.from_pandas_edgelist(edge_list, source = 'SIRE', target = 'ID',create_using=nx.DiGraph()) #load in nx
    
    founders = ped_sub_sort.loc[ped_sub_sort['sire'] == '0']
    for founder in founders.iloc[:,0]:
      des = nx.descendants(G, source=founder)
      ped_sub_sub = PyAGH.selectPed(data=ped_sub_sort,id=list(map(int, des)),generation=1)

    #Calculate additive relationship matrix
      A = PyAGH.makeA(ped_sub_sub,Sparse=False)
    
    #Loop through matrix and take out pairs
      N = len(A[1]) #Number of individuals
    
      for k in range(N):
        for l in range(k+1, N):
          if A[0][k,l] != 0:
            f.write(str(A[1][k]) + ' ' + str(A[1][l]) + ' ' + str(A[0][k,l]) + '\n')


f.close()

