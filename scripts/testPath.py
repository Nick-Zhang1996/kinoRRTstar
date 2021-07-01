import numpy as np
import matplotlib.pyplot as plt
from math import sin,cos

def getTrajectory(t, der=0):
    T = 10
    R = 1.0
    A = 0.2
    e = 0.01
    x = lambda t: np.cos(t/T * 2*np.pi) * (R + A*np.sin(t/T*4* 2*np.pi))
    y = lambda t: np.sin(t/T * 2*np.pi) * (R + A*np.sin(t/T*4* 2*np.pi))
    z = lambda t: -0.3

    if (der == 0):
        return np.array((x(t),y(t),z(t)))

    deri = lambda fun,t: (fun(t+e) - fun(t-e))/(2*e)
    if (der == 1):
        return np.array((deri(x,t),deri(y,t),deri(z,t)))

    dderi = lambda fun,t: (fun(t+e) -2*fun(t) + fun(t-e))/(e*e)
    if (der == 2):
        return np.array((dderi(x,t),dderi(y,t),dderi(z,t)))

def plotDir(loc,direction):
    loc = np.array(loc)
    direction = np.array(direction)
    p0 = loc
    p1 = loc + direction
    plt.plot( [p0[0],p0[0] + p1[0]], [p1[1],p0[1]+p1[1]] ,'*')
    


tt = np.linspace(0,10,1000)
xyz_vec = []
dxyz_vec = []
ddxyz_vec = []
for t in tt:
    xyz_vec.append(getTrajectory(t))
    dxyz_vec.append(getTrajectory(t,der=1))
    ddxyz_vec.append(getTrajectory(t,der=2))

xyz_vec = np.array(xyz_vec)
dxyz_vec = np.array(dxyz_vec)
ddxyz_vec = np.array(ddxyz_vec)

plt.plot(xyz_vec[:,0], xyz_vec[:,1])

index = 20
loc = (xyz_vec[index,0], xyz_vec[index,1])
dire = (dxyz_vec[index,0], dxyz_vec[index,1])
plotDir(loc, dire)
'''
dire = (ddxyz_vec[index,0], ddxyz_vec[index,1])
plotDir(loc, dire)
'''
plt.show()



