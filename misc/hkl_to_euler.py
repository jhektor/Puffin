import numpy as np
import sys

# Plane normal
h = int(sys.argv[1])
k = int(sys.argv[2])
l = int(sys.argv[3])
# Direction
u = int(sys.argv[4])
v = int(sys.argv[5])
w = int(sys.argv[6])

n = 1/(np.sqrt(h**2+k**2+l**2))*np.array([h,k,l])
b = 1/(np.sqrt(u**2+v**2+w**2))*np.array([u,v,w])
tn = np.cross(n,b)
t = tn/np.linalg.norm(tn)
U = np.stack((b,t,n),axis=1)

# From matrix to Euler angles
# special case if PHI=0
if np.abs(U[2,2]-1)<1e-5:
    PHI=0
    phi1 = np.arctan2(U[0,1],U[0,0])/2.
    phi2 = -phi1
else:
    PHI = np.arccos(U[2,2])
    phi2 = np.arctan2(U[0,2]/np.sin(PHI),U[1,2]/np.sin(PHI))
    phi1 = np.arctan2(U[2,0]/np.sin(PHI),-U[2,1]/np.sin(PHI))



# From radians to degrees
PHI = np.rad2deg(PHI)
phi1 = np.rad2deg(phi1)
phi2 = np.rad2deg(phi2)

print 'Euler angles corresponding to (',h,k,l,')[',u,v,w,']:', phi1,PHI,phi2
