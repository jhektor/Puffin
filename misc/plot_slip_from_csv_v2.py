import numpy as np
import glob
import sys

import matplotlib.pyplot as plt

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


#Scalar postprocessors
scalar_csv = sys.argv[-1]# slip_csv[1][:-19]+'.csv'
pps = np.genfromtxt(scalar_csv,delimiter=',')
time = pps[1:,0]
strain_zz = pps[1:,-2]


plt.figure()
for j,col in enumerate(pps[1:,1:-2].T):
    plt.plot(strain_zz,col,label=j)
    print 'slip system:', j+1, ' max slip: ', np.max(np.abs(col))
plt.legend()
plt.show()
# r = np.array(res_csv)
# np.savetxt(slip_csv[1][:-19]+'_res.csv',r.T,delimiter=',')
