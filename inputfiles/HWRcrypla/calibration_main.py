from subprocess import call #interface to command LangevinNoisePositive
import scipy.optimize as spo
import scipy.io as sio
import numpy as np

import matplotlib.pyplot as plt


plt.ion() #interactive plotting

fig = plt.figure()
ax1 = fig.add_subplot(111)
# ax1.set_title('001 loading')
# ax1.set_ylim([0,35])
# ax2 = fig.add_subplot(122)
# ax2.set_title('100 loading')
# ax2.set_ylim([0,35])

#Define the objective function, minimize the 2-norm or simulation-experiment
def calibfcn(x, pltstring = '--'):
    error = 0

    G0 = [8.5, 10.4, 5.1, 5.1, 8.5, 10.4, 5.1, 5.1, 4.3, 4.5, 4.5, 5.6, 4.3, 4.5, 4.5, 5.6, 7.4, 7.4, 15, 15, 6.6, 6.6, 6.6, 6.6, 12, 12, 12, 12, 12, 12, 12, 12]
    for loadcase in [0]:
        #Read experimental data from mat files
        if loadcase is 0:
            data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/Kariya010.mat'
            eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=90 Materials/elasticity_tensor/euler_angle_3=0' #changes the rotation of crystal to 010 along loading in input file
            G0[2] = 2#5.1
            G0[3] = 2#5.1
            G0[5] = 5#10.4
            G0[6] = 2#5.1
            G0[7] = 2#5.1
            G0[15] = 3#5.6
            G0[16] =4# 7.4
            G0[24] = 6#12
            G0[25] = 6#12

            # data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/5compDirTest_001.mat'
            # # data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/philippi2016extra001.mat'
            # eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=0 Materials/elasticity_tensor/euler_angle_3=0' #changes the rotation of crystal to 001 along loading in input file
            # # eul2 = 'Materials/elasticity_tensor/euler_angle_1=60 Materials/elasticity_tensor/euler_angle_2=90' #changes the rotation of crystal to 001 along loading in input file
        elif loadcase is 1:
            data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/Kariya110.mat'
            eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=90 Materials/elasticity_tensor/euler_angle_3=45' #changes the rotation of crystal to 010 along loading in input file            # data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/5compDirTest_100.mat'
            # # data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/philippi2016extra110v2.mat'
            # eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=90 Materials/elasticity_tensor/euler_angle_3=90' #changes the rotation of crystal to 100 along loading in input file
            # # eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=90 Materials/elasticity_tensor/euler_angle_3=45' #changes the rotation of crystal to 100 along loading in input file
            # # eul2 = 'Materials/elasticity_tensor/euler_angle_1=30 Materials/elasticity_tensor/euler_angle_2=60 Materials/elasticity_tensor/euler_angle_3=75' #
        data = sio.loadmat(data_file)
        strain_exp = data['xx'][:,0]
        stress_exp = data['yy'][:,0]#*1e-6 #in MPa

        x0 = [50, 20, 0.0001, 6] # B, Q, gamma0, m
        xs = x*x0

        #Set up moose input file
        inputfile = '1element_HWR.i '
        # names of properties to calibrate
        state_var_rate_B_name = 'UserObjects/state_var_evol_rate/B='
        slip_res_Q_name = 'UserObjects/slip_resistance/Q='
        gamma0_name = 'UserObjects/slip_rate/gamma0='
        m_name = 'UserObjects/slip_rate/m='
        G0_name = 'UserObjects/slip_resistance/G0='

        B_str = "'" + str(xs[0]) + "'"
        Q_str = str(xs[1])
        gamma0_str = str(xs[2])
        m_str = str(xs[3])
        G0_str = str(G0)
        G0_str = G0_str.replace('[',"'")
        G0_str = G0_str.replace(']',"'")
        G0_str = G0_str.replace(', '," ")
        #Run moose simulation
        print 'Load case:', loadcase
        print "\033[95mCurrent material parameters:" + "\033[95m{}\033[0m".format(xs)
        runcmd = 'mpirun -n 1 ../../puffin-opt -i ' + inputfile + state_var_rate_B_name + B_str + ' ' + slip_res_Q_name + Q_str + ' ' + gamma0_name + gamma0_str + ' ' + m_name + m_str + ' ' + eul2 +' ' + G0_name + G0_str + ' > mooselog.txt'
        print 'Running this command:\n' + runcmd + "\033[0m"
        call(runcmd, shell=True)

        #Get stress strain curve from csv file
        # aa = np.recfromcsv('calibrationSn.csv')
        aa = np.loadtxt('hwr.csv',delimiter = ',', skiprows = 1)
        # idx = (np.abs(-aa[:,-3] - 0.12)).argmin()
        #idx = -1
        strain_sim = -aa[1:,-2] #eps_zz
        stress_sim = -aa[1:,-1] #sigma_zz in MPa (compression positive)
        print error
        if np.max(strain_sim) < 0.048: #this means the simulation failed ???
            error += 2000
        else:
            #Interpolate experimental values to simulated times
            stress_exp_interp = np.interp(strain_sim,strain_exp,stress_exp)
            #Calculate error
            # error += np.linalg.norm((stress_sim-stress_exp_interp)/stress_exp_interp)
            error += np.linalg.norm((stress_exp_interp-stress_sim)/stress_sim)

        if loadcase is 0:
            # error = np.linalg.norm((stress_sim-stress_exp_interp)/stress_exp_interp)
            pltstring = '-'
            ax1.plot(strain_exp,stress_exp,'ko')
            ax1.plot(strain_sim,stress_sim,pltstring)
        elif loadcase is 1:
            pltstring = '--'
            ax1.plot(strain_exp,stress_exp,'go')
            ax1.plot(strain_sim,stress_sim,pltstring)

        plt.pause(0.05)
        # plt.pause(5)

    print "\033[91mError is: \033[00m"+"\033[91m {}\033[00m".format(error)

    return error
# Minimize the objective function
# x = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
# bounds = ((0,None),(0,None),(0,None),(0,None),(0,None))
bounds = ((0, None),(0,None), (0,None),(0,20))

# mpar = [0.001, 6, 7e-3, 8]
# mpar = 67*[1]
mpar = [1, 1, 1,1] #B, Q, gamma0, m

results = spo.minimize(calibfcn,mpar,bounds=bounds)
print mpar
if not results.success:
    print results.message
else:
    print "Successful optimization!, %5d, iterations" % (results.nit)
#Run simulation with the calibrated parameters
calibfcn(results.x,pltstring='-')
# calibfcn(mpar,pltstring='-')
# plt.pause()
# ax.plot(strain_exp,stress_exp,strain_sim,stress_sim,strain_sim,stress_exp_interp)

plt.show(block=True)
# calibfcn(mpar,(strain_exp, stress_exp))
