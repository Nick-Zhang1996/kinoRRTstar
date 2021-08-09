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
from window import createWindow


# ------  Load data ------
logFilename = "./log.p"
output = open(logFilename,'rb')
data = pickle.load(output)
output.close()

logFilename = "./rrt_traj.p"
output = open(logFilename,'rb')
waypoints = pickle.load(output)
output.close()

data = np.array(data)
print("log length %.1f seconds"%(data.shape[0]/119.88))

skip = 0
t = data[skip:,0] - data[0,0]
x = data[skip:,1]
y = data[skip:,2]
z = data[skip:,3]
rx = data[skip:,4]
ry = data[skip:,5]
rz = data[skip:,6]


# ------ Misc Analytics ------
print_ok("desired")
desired_max_speed = (np.max(waypoints[:,4]**2 + waypoints[:,5]**2 + waypoints[:,6]**2))**0.5
desired_max_acc = (np.max(waypoints[:,7]**2 + waypoints[:,8]**2 + waypoints[:,9]**2))**0.5
print_info("trajectory time %.2f" %(waypoints[-1,0]-waypoints[0,0]))
print_info("max speed : %.1f m/s "%(desired_max_speed))
print_info("max acc : %.1f m/s "%(desired_max_acc))

print_ok("actual:")
dx = np.diff(x)/(t[1]-t[0])
dy = np.diff(y)/(t[1]-t[0])
ddx = np.diff(dx)/(t[1]-t[0])
ddy = np.diff(dy)/(t[1]-t[0])
max_speed = np.max((dx*dx + dy*dy))**0.5
max_acc = np.max((ddx*ddx + ddy*ddy))**0.5
print_info("trajectory time %.2f" %(t[-1]))
print_info("(?)max speed : %.1f m/s "%(max_speed))
print_info("(?)max acc : %.1f m/s "%(max_acc))
'''
print_info("speed")
plt.plot((dx*dx + dy*dy)**0.5)
plt.show()
print_info("acc")
plt.plot((ddx*ddx + ddy*ddy)**0.5,'o')
plt.show()
print_info("roll/pitch")
plt.plot(rx*180.0/np.pi)
plt.plot(ry*180.0/np.pi)
plt.show()
'''

start_node = Node(-1.8, 0.6, -0.6)
goal_node = Node(1.7,0,-2)

world = World(-2,2,-0.7,1,-2.5,0)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)
obstacles = []
obstacles += createWindow((-1,0,0), 0.5, True, world)
obstacles += createWindow((0.5,-1,0), 0.5, False, world)
# ------  Plot World ------
for obstacle in obstacles:
    visual.addObstacle(obstacle)
ax = visual.visualizeWorld(show=False)
# ------ plot desired trajectory ( waypoints) ------
ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'b', label='desired')

# ------ plot actual trajectory
ax.plot(-x, y, -z, color='r', label='actual')
ax.set_xlabel('-x')
ax.set_ylabel('-y')
ax.set_zlabel('altitude (-z)')
ax.legend()
plt.show()


