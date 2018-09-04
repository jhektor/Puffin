from subprocess import call #interface to command LangevinNoisePositive
import scipy.optimize as spo
import scipy.io as sio
import numpy as np

import matplotlib.pyplot as plt


plt.ion() #interactive plotting

fig = plt.figure()
ax1 = fig.add_subplot(131)
ax1.set_title('010 loading')
ax1.set_ylim([0,35])
ax2 = fig.add_subplot(132)
ax2.set_title('100 loading')
ax2.set_ylim([0,35])
ax3 = fig.add_subplot(133)
ax3.set_title('101 loading')
ax3.set_ylim([0,35])

#Define the objective function, minimize the 2-norm or simulation-experiment
def calibfcn(x, pltstring = '--'):
    error = 0
    
    for loadcase in [2]:
        #Read experimental data from mat files
        if loadcase is 0:
            data_file = '/home/viktor/projects/Puffin/inputfiles/calibration_crypla/data/5compDirTest_010.mat'
            eul2 = 'Materials/elasticity_tensor/euler_angle_1=90 Materials/elasticity_tensor/euler_angle_2=90 Materials/elasticity_tensor/euler_angle_3=0' #changes the rotation of crystal to 001 along loading in input file            
        elif loadcase is 1:
            data_file = '/home/viktor/projects/Puffin/inputfiles/calibration_crypla/data/5compDirTest_100.mat'
            eul2 = 'Materials/elasticity_tensor/euler_angle_1=180 Materials/elasticity_tensor/euler_angle_2=90 Materials/elasticity_tensor/euler_angle_3=90' #changes the rotation of crystal to 100 along loading in input file
        elif loadcase is 2:
            data_file = '/home/viktor/projects/Puffin/inputfiles/calibration_crypla/data/5compDirTest_101.mat'
            eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=45 Materials/elasticity_tensor/euler_angle_3=90' #changes the rotation of crystal to 101 along loading in input file
        data = sio.loadmat(data_file)
        strain_exp = data['xx'][:,0]
        stress_exp = data['yy'][:,0]*1e-6 #in MPa

        #Set up moose input file
        inputfile = '1element_calib_removed.i '
        # names of properties to calibrate
        slip_rate_props_name = 'UserObjects/slip_rate_gss/flowprops='
        state_var_props_name = 'UserObjects/state_var_gss/group_values='
        state_var_rate_h0_name = 'UserObjects/state_var_evol_rate_comp_gss/h0_group_values='
        state_var_rate_tauSat_name = 'UserObjects/state_var_evol_rate_comp_gss/tauSat_group_values='
        state_var_rate_hardeningExponent_name = 'UserObjects/state_var_evol_rate_comp_gss/hardeningExponent_group_values='

        # initial values (from Darbandi2012)
        slip_rate_props_vals = [1, 20, 0.001, 0.05] #start_ss end_ss gamma0 1/m m = 20??
        state_var_props_vals = [x[0]]*10# initial slip resistance values of each ss
        state_var_rate_h0_vals =  [x[1]]*10 #  h0 of each ss
        state_var_rate_tauSat_vals = [x[2]]*10 # tau saturation of each ss
        state_var_rate_hardeningExponent_vals = [x[3]]*10 # the hardening exponent c

        slip_rate_props = '\''+" ".join(str(x) for x in slip_rate_props_vals)+'\' '
        state_var_props = '\''+" ".join(str(x) for x in state_var_props_vals)+'\' '
        state_var_rate_h0 = '\''+" ".join(str(x) for x in state_var_rate_h0_vals)+'\' '
        state_var_rate_tauSat = '\''+" ".join(str(x) for x in state_var_rate_tauSat_vals)+'\' '
        state_var_rate_hardeningExponent = '\''+" ".join(str(x) for x in state_var_rate_hardeningExponent_vals)+'\' '


        #Run moose simulation
        print 'Load case:', loadcase
        print "\033[94mCurrent material parameters [MPa]:" + "\033[94m{}\033[0m".format(x*160.217662)
        print "\033[95mCurrent material parameters:" + "\033[95m{}\033[0m".format(x)
        runcmd = 'mpirun -n 1 ../../puffin-opt -i ' + inputfile + slip_rate_props_name + slip_rate_props + state_var_props_name + state_var_props + state_var_rate_h0_name + state_var_rate_h0 + state_var_rate_tauSat_name + state_var_rate_tauSat + state_var_rate_hardeningExponent_name + state_var_rate_hardeningExponent + eul2 + ' > mooselog.txt'
        print 'Running this command:\n' + runcmd + "\033[0m"
        call(runcmd, shell=True)

        #Get stress strain curve from csv file
        # aa = np.recfromcsv('calibrationSn.csv')
        aa = np.loadtxt('calibrationSn.csv',delimiter = ',', skiprows = 1)
        # idx = (np.abs(-aa[:,-3] - 0.12)).argmin()
        #idx = -1
        strain_sim = -aa[:,-3] #eps_yy
        stress_sim = -aa[:,-1]*160.217662 #sigma_yy in MPa (compression positive)

        if np.max(strain_sim) < 0.048: #this means the simulation failed ???
            error += 20
        else:
            #Interpolate experimental values to simulated times
            stress_exp_interp = np.interp(strain_sim,strain_exp,stress_exp)
            #Calculate error
            error += np.linalg.norm((stress_sim-stress_exp_interp)/stress_exp_interp)

        if loadcase is 0:
            # error = np.linalg.norm((stress_sim-stress_exp_interp)/stress_exp_interp)
            ax1.plot(strain_exp,stress_exp,'ko')
            ax1.plot(strain_sim,stress_sim,pltstring)
        elif loadcase is 1:
            ax2.plot(strain_exp,stress_exp,'ko')
            ax2.plot(strain_sim,stress_sim,pltstring)
        elif loadcase is 2:
            ax3.plot(strain_exp,stress_exp,'ko')
            ax3.plot(strain_sim,stress_sim,pltstring)

        plt.pause(0.05)

    print "\033[91mError is: \033[00m"+"\033[91m {}\033[00m".format(error)

    return error
# Minimize the objective function
# x = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
bounds = ((0, None),)*4    
#bounds = ((1e-3, 1e-1),) + ((5e-3, 1e-1),) + ((1e-2, 0.8),) + ((0.15, 0.75),) + ((1e-1, 0.75),) + ((1.5, 3),) + (((1e-3, 1e-1),) + ((5e-3, 1e-1),))*2 # set bounds for all variables

#Initial values of parameter to calibrate (from Darbandi2012) [h0, ss, s0] used to scale x
#Add more s0 varibles. Take initial values from table 8.4 and 8.5 in Darbandis thesis
# mpar = [0.468, 0.075, 0.053, 0.0268, 0.0649, 0.0281, 0.0350, 0.0318, 0.0462, 0.0936, 0.0412, 0.0749]
# mpar = [0.44401723, 0.01013415, 0.05331709, 0.02615524, 0.06898327, 0.02759656, 0.03323911, 0.01238013, 0.04818832, 0.09332103, 0.6712851, 0.02919173]
# mpar = [0.44401723, 0.075, 0.05331709, 0.02615524, 0.064912, 0.0281, 0.034952, 0.031832, 0.046187, 0.093623, 0.041194, 0.074898]
# x = mpar*x

mpar = np.array([0.14, 0.8, 0.24, 2])

#after 101 opt.
#[0.00923856 0.19937991 0.07224225 2.54105146]

results = spo.minimize(calibfcn,mpar,bounds=bounds)
print mpar*160.217662 #result
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
