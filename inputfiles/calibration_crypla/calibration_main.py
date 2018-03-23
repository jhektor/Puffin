from subprocess import call #interface to command LangevinNoisePositive
import scipy.optimize as spo
import scipy.io as sio
import numpy as np

import matplotlib.pyplot as plt


plt.ion() #interactive plotting

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax1.set_title('001 loading')
ax1.set_ylim([0,30])
ax2 = fig.add_subplot(122)
ax2.set_title('110 loading')
ax2.set_ylim([0,30])

#Define the objective function, minimize the 2-norm or simulation-experiment
def calibfcn(x, finished = False):
    error = 0
    for loadcase in [0,1]:
        #Read experimental data from mat files
        if loadcase is 0:
            data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/philippi2016extra001.mat'
            eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=0 Materials/elasticity_tensor/euler_angle_3=0' #changes the rotation of crystal to 001 along loading in input file
            # eul2 = 'Materials/elasticity_tensor/euler_angle_1=60 Materials/elasticity_tensor/euler_angle_2=90' #changes the rotation of crystal to 001 along loading in input file
        elif loadcase is 1:
            data_file = '/home/johan/projects/Puffin/inputfiles/calibration_crypla/data/philippi2016extra110v2.mat'
            eul2 = 'Materials/elasticity_tensor/euler_angle_1=0 Materials/elasticity_tensor/euler_angle_2=90 Materials/elasticity_tensor/euler_angle_3=45' #changes the rotation of crystal to 110 along loading in input file
            # eul2 = 'Materials/elasticity_tensor/euler_angle_1=45 Materials/elasticity_tensor/euler_angle_2=60' #changes the rotation of crystal to 110 along loading in input file
        data = sio.loadmat(data_file)
        strain_exp = data['xx'][:,0]
        stress_exp = data['yy'][:,0]*1e-6 #in MPa

        #Set up moose input file
        inputfile = '1element_calib.i '
        # names of properties to calibrate
        slip_rate_props_name = 'UserObjects/slip_rate_gss/flowprops='
        state_var_props_name = 'UserObjects/state_var_gss/group_values='
        state_var_rate_props_name = 'UserObjects/state_var_evol_rate_comp_gss/hprops='

        # initial values (from Darbandi2012)
        slip_rate_props_vals = [1, 32, 0.001, 0.05] #start_ss end_ss gamma0 1/m
        state_var_rate_props_vals = [1.4, x[0], x[1], 2] #qab h0 ss c see eq (9) in Zhao 2017
        state_var_vals = [x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11]] #slip resistance

        slip_rate_props = '\''+" ".join(str(x) for x in slip_rate_props_vals)+'\' '
        state_var_props = '\''+" ".join(str(x) for x in state_var_vals)+'\' '
        state_var_rate_props = '\''+" ".join(str(x) for x in state_var_rate_props_vals)+'\' '

        #Run moose simulation
        print 'Load case:', loadcase
        print "\033[95mCurrent material parameters:" + "\033[95m{}\033[0m".format(x)
        runcmd = 'mpirun -n 1 ../../puffin-opt -i ' + inputfile + slip_rate_props_name + slip_rate_props + state_var_props_name + state_var_props + state_var_rate_props_name + state_var_rate_props + eul2 + ' > mooselog.txt'
        print 'Running this command:\n' + runcmd + "\033[0m"
        call(runcmd, shell=True)

        #Get stress strain curve from csv file
        # aa = np.recfromcsv('calibrationSn.csv')
        aa = np.loadtxt('calibrationSn.csv',delimiter = ',', skiprows = 1)
        strain_sim = -aa[:,2] #eps_yy
        stress_sim = -aa[:,-1]*160.217662 #sigma_yy in MPa (compression positive)

        if np.max(strain_sim) < 0.12: #this means the simulation failed
            error += 20
        else:
            #Interpolate experimental values to simulated times
            stress_exp_interp = np.interp(strain_sim,strain_exp,stress_exp)
            #Calculate error
            error += np.linalg.norm((stress_sim-stress_exp_interp)/stress_exp_interp)

        if loadcase is 0:
            # error = np.linalg.norm((stress_sim-stress_exp_interp)/stress_exp_interp)
            ax1.plot(strain_exp,stress_exp,'ko')
            ax1.plot(strain_sim,stress_sim)
        elif loadcase is 1:
            ax2.plot(strain_exp,stress_exp,'ko')
            ax2.plot(strain_sim,stress_sim)

        plt.pause(0.05)

    print "\033[91mError is: \033[00m"+"\033[91m {}\033[00m".format(error)

    return error
# Minimize the objective function
# x = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
bounds = ((5e-3, None),(5e-3,None),(5e-3,None),(5e-3,None),(5e-3,None),(5e-3,None),(5e-3,None),(5e-3,None),(5e-3,None),(5e-3,None),(5e-3,None),(5e-3,None))

#Initial values of parameter to calibrate (from Darbandi2012) [h0, ss, s0] used to scale x
#Add more s0 varibles. Take initial values from table 8.4 and 8.5 in Darbandis thesis
# mpar = [0.468, 0.075, 0.053, 0.0268, 0.0649, 0.0281, 0.0350, 0.0318, 0.0462, 0.0936, 0.0412, 0.0749]
# mpar = [0.44401723, 0.01013415, 0.05331709, 0.02615524, 0.06898327, 0.02759656, 0.03323911, 0.01238013, 0.04818832, 0.09332103, 0.6712851, 0.02919173]
mpar = [0.44401723, 0.075, 0.05331709, 0.02615524, 0.064912, 0.0281, 0.034952, 0.031832, 0.046187, 0.093623, 0.041194, 0.074898]
# x = mpar*x
result = spo.minimize(calibfcn,mpar,bounds=bounds)
print result

#Run simulation with the calibrated parameters
# calibfcn(results.x)
# calibfcn(mpar)
# plt.pause()
# ax.plot(strain_exp,stress_exp,strain_sim,stress_sim,strain_sim,stress_exp_interp)




plt.show(block=True)
# calibfcn(mpar,(strain_exp, stress_exp))
