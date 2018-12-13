import numpy as np
import glob
import sys

import matplotlib.pyplot as plt

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx



#Find slip, eta1 eta2 csv files first file is empty for unknown reason...
slip_csv = sorted(glob.glob('*slip_line*.csv'))[1:]
eta1_csv = sorted(glob.glob('*eta1_line*.csv'))[1:]
eta2_csv = sorted(glob.glob('*eta2_line*.csv'))[1:]

#Scalar postprocessors
scalar_csv = slip_csv[1][:-19]+'.csv'
pps = np.genfromtxt(scalar_csv,delimiter=',')
time = pps[2:,0]
plot_times_rec = [25,50,100]
plot_times_sim = []
idxs = []
for pt in plot_times_rec:
    tt,idx = find_nearest(time,pt)
    plot_times_sim.append(tt)
    idxs.append(idx)

print plot_times_sim
print idxs

# x-coordinates (always the same)
x = np.genfromtxt(slip_csv[1],delimiter=',')[:,0]

res_csv = [x]

plt.figure()
for j,i in enumerate(idxs):
    e1 = eta1_csv[i]
    e2 = eta2_csv[i]
    s = slip_csv[i]
    eta1 = np.genfromtxt(e1,delimiter=',')[:,4]
    eta2 = np.genfromtxt(e2,delimiter=',')[:,4]
    slip = np.genfromtxt(s,delimiter=',')[:,4]
    result = (eta1**2+eta2**2)*slip
    result[result<0] = 0
    res_csv.append(result)
    plt.plot(x,result,label=plot_times_rec[j])

plt.legend()
r = np.array(res_csv)
np.savetxt(slip_csv[1][:-19]+'_res.csv',r.T,delimiter=',')
