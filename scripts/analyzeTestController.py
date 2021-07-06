import pickle
import matplotlib.pyplot as plt
import numpy as np
from common import *
from time  import time,sleep
import os
import sys
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src/')
sys.path.append(base_dir)
from kinoRRT import *
from WorldVisualization import WorldVisualization

logFilename = "./log.p"
output = open(logFilename,'rb')
data = pickle.load(output)
output.close()

data = np.array(data)
print("log length %.1f seconds"%(data.shape[0]/119.88))

t = data[:,0]
x = data[:,1]
y = data[:,2]
z = data[:,3]
rx = data[:,4]
ry = data[:,5]
rz = data[:,6]


# plot desired trajectory
T = 10
t = np.linspace(0,5)
R = 1.0
A = 0.2
x_t = lambda t: np.cos(t/T * 2*np.pi) * (R + A*np.cos(t/T*4* 2*np.pi))
y_t = lambda t: np.sin(t/T * 2*np.pi) * (R + A*np.cos(t/T*4* 2*np.pi))
z_t = lambda t: -0.5

x_des = []
y_des = []
z_des = []
for t in np.linspace(0,10):
    x_des.append(x_t(t))
    y_des.append(y_t(t))
    z_des.append(z_t(t))

x_des = np.array(x_des)
y_des = np.array(y_des)
z_des = np.array(z_des)

print("yaw")
plt.plot(rz*180.0/np.pi)
plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot([0,0],[0,0],[-1,1],'white')
ax.plot(-x_des, y_des, -z_des,'b--',label='desired')
ax.set_xlabel("-x")
ax.set_ylabel("y")
ax.set_zlabel("-z")




# Plot 3d.
'''
fig = plt.figure()
ax = fig.gca(projection='3d')
'''

ax.plot(-x, y, -z, color='r', label='actual')
ax.set_xlabel('-x')
ax.set_ylabel('-y')
ax.set_zlabel('altitude (-z)')
ax.legend()
plt.show()
