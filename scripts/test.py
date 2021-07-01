import numpy as np
from math import sin,cos,radians,degrees
from scipy.spatial.transform import Rotation as R

def getPassiveR(r,p,y):
    R = [[cos(p) * cos(y), cos(p)*sin(y), -sin(p)],[sin(r)*sin(p)*cos(y) - cos(r)*sin(y), sin(r)*sin(p)*sin(y) + cos(r)*cos(y), cos(p)*sin(r)], [cos(r)*sin(p)*cos(y)+sin(r)*sin(y), cos(r)*sin(p)*sin(y)-sin(r)*cos(y), cos(p)*cos(r)]]
    R = np.array(R)
    return R

def getAciveR(r,p,y):
    return getPassiveR(r,p,y).T


#z = np.array([0,0,1])
#print(getAciveR(radians(90),0,radians(90)) @ z)
Fdes = np.array([0,1,0])
T = np.array([0,0,-1])

axis = ( Fdes/np.linalg.norm(Fdes) + T / np.linalg.norm(T))
r_vec = axis / np.linalg.norm(axis) * np.pi
r = R.from_rotvec(r_vec)

print("error")
print(r.as_matrix() @ T - Fdes)
print(r.as_euler('ZYX')*180/np.pi)

