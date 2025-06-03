#!/usr/bin/python

import pandas as pd

#The biggest family
res = pd.read_table('/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/relatedness_coeffs_biggest_IMID.txt.qz', sep = ' ', header = None)

res.columns = ['id1', 'id2', 'r']

res = res[~res[['id1', 'id2']].apply(frozenset, axis=1).duplicated()]

export_csv = res.to_csv(r'/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/relatedness_coeffs_biggest_no_duplicates_IMID.csv.gz', compression = 'gzip', sep = ' ', index = None, header = True)


#The rest of the families
res = pd.read_table('/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/relatedness_coeffs_less_IMID.txt', sep = ' ', header = None)

res.columns = ['id1', 'id2', 'r']

res = res[~res[['id1', 'id2']].apply(frozenset, axis=1).duplicated()]

export_csv = res.to_csv(r'/ngc/projects2/predict_r/research/projects/0015_Genealogy_across_diseases/Generated_Data/relatedness_coeffs_less_no_duplicates_IMID.csv.gz', compression = 'gzip', sep = ' ', index = None, header = True)
