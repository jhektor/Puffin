from subprocess import call #interface to command LangevinNoisePositive
import scipy.optimize as spo
import scipy.io as sio
import numpy as np

import matplotlib.pyplot as plt


plt.ion() #interactive plotting

fig = plt.figure()
ax1 = fig.add_subplot(111)
# ax1.set_title('001 loading')
ax1.set_ylim([0,35])
# ax2 = fig.add_subplot(122)
# ax2.set_title('100 loading')
# ax2.set_ylim([0,35])

#Define the objective function, minimize the 2-norm or simulation-experiment
def calibfcn(x, pltstring = '--'):
    error = 0

    for loadcase in [0]:
        #Read experimental data from mat files
        if loadcase is 0:
            # data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/5compDirTest_001.mat'
            data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/philippi2016extra001.mat'
            eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=0 Materials/elasticity_tensor/euler_angle_3=0' #changes the rotation of crystal to 001 along loading in input file
            # eul2 = 'Materials/elasticity_tensor/euler_angle_1=60 Materials/elasticity_tensor/euler_angle_2=90' #changes the rotation of crystal to 001 along loading in input file
        elif loadcase is 1:
            # data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/5compDirTest_100.mat'
            data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/philippi2016extra110v2.mat'
            # eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=90 Materials/elasticity_tensor/euler_angle_3=-45' #changes the rotation of crystal to 100 along loading in input file
            eul2 = 'Materials/elasticity_tensor/euler_angle_1=90 Materials/elasticity_tensor/euler_angle_2=90 Materials/elasticity_tensor/euler_angle_3=45' #changes the rotation of crystal to 100 along loading in input file
            # eul2 = 'Materials/elasticity_tensor/euler_angle_1=30 Materials/elasticity_tensor/euler_angle_2=60 Materials/elasticity_tensor/euler_angle_3=75' #
        data = sio.loadmat(data_file)
        strain_exp = data['xx'][:,0]
        stress_exp = data['yy'][:,0]*1e-6 #in MPa

        #Set up moose input file
        inputfile = '1element_HWR.i '
        # names of properties to calibrate
        slip_rate_props_gamma0_name = 'UserObjects/slip_rate/gamma0='
        slip_rate_props_m_name = 'UserObjects/slip_rate/m='
        state_var_props_name = 'UserObjects/state_var/group_values='
        state_var_rate_B_name = 'UserObjects/state_var_evol_rate/B='
        slip_res_Q_name = 'UserObjects/slip_resistance/Q='

        # initial values (from Darbandi2012)
        scaling = [0.001, 6, 20]
        scaling.extend([8 for i in range(32)])
        scaling.extend([7e-3 for i in range(32)])
        slip_rate_props_gamma0 = x[0]*scaling[0] #0.001
        slip_rate_props_m = x[1]*scaling[1] #6
        slip_res_Q = x[2]*scaling[2] #20
        state_var_rate_B =  x[3:35]*scaling[3:35] #8.0
        state_var_props = x[35:]*scaling[35:] #7e-3 # initial slip resistance values of each ss

        state_var_str = '\''+" ".join(str(sv) for sv in state_var_props)+'\' '
        state_var_rate_B_str = '\''+" ".join(str(sv) for sv in state_var_rate_B)+'\' '


        #Run moose simulation
        print 'Load case:', loadcase
        print "\033[95mCurrent scaling parameters:" + "\033[95m{}\033[0m".format(x)
        print "\033[95mCurrent material parameters:" + "\033[95m{}\033[0m".format(x*scaling)
        runcmd = 'mpirun -n 1 ../../puffin-opt -i ' + inputfile + slip_rate_props_gamma0_name + str(slip_rate_props_gamma0)+ ' ' + slip_rate_props_m_name + str(slip_rate_props_m)+ ' ' + state_var_props_name + state_var_str+ ' ' + state_var_rate_B_name + state_var_rate_B_str + ' ' + slip_res_Q_name + str(slip_res_Q) + ' ' + eul2 + ' > mooselog.txt'
        print 'Running this command:\n' + runcmd + "\033[0m"
        call(runcmd, shell=True)

        #Get stress strain curve from csv file
        # aa = np.recfromcsv('calibrationSn.csv')
        aa = np.loadtxt('hwr.csv',delimiter = ',', skiprows = 1)
        # idx = (np.abs(-aa[:,-3] - 0.12)).argmin()
        #idx = -1
        strain_sim = -aa[:,-3] #eps_zz
        stress_sim = -aa[:,-1] #sigma_zz in MPa (compression positive)

        if np.max(strain_sim) < 0.048: #this means the simulation failed ???
            error += 20
        else:
            #Interpolate experimental values to simulated times
            stress_exp_interp = np.interp(strain_sim,strain_exp,stress_exp)
            #Calculate error
            error += np.linalg.norm((stress_sim-stress_exp_interp)/stress_exp_interp)

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
bounds = ((0, None),)*67

# mpar = [0.001, 6, 7e-3, 8]
# mpar = 67*[1]
mpar = np.random.rand(67)

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
