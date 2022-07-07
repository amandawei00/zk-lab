import numpy as np
import pandas as pd
import csv

# import experimental data
data = pd.read_csv('../data/redx2009_full.csv', delimiter='\t', header=0, comment='#')

# import farid's solutions
farid = pd.read_csv('mve_dis_farid.csv', header=0, delimiter='\t')

print(data)
print(farid)

exp  = np.array(data['redx'])
err  = np.array(data['err'])
the  = np.array(farid['redx'])

err  = 0.01 * exp * err
chi2 = np.sum(np.power((the - exp)/err, 2))/len(the)
print(chi2)  
