
import os
import sys
import numpy as np
from time import time
import matplotlib.pyplot as plt

base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../build/')
sys.path.append(base_dir)
from kinoRRT import *

from WorldVisualization import WorldVisualization
from common import *
from window import createWindow


# Original big test world
world = World(-10,10,-10,10,-10,10)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)

obs1 = Box(-4,-2, -10,0, -10,10)
obs2 = Box(-4,-2, 0,10, -10,0)
obs3 = Box(2,4, 0,10, -10,10)
obs4 = Box(2,4, -10,0, 0,10)
start_node = Node(-8, 2, -3)
goal_node = Node(8, -2, 3)
obstacles = [obs1, obs2, obs3, obs4]

for obstacle in obstacles:
    world.addObstacle(obstacle)
    visual.addObstacle(obstacle)

visual.visualizeWorld(show=True)
duration = 30.0
print("solving with SST, time limit = %.2f sec"%(duration))
sst = SST(world, start_node, goal_node, duration)
#sst.solve()
count = sst.getWaypointCount()
waypoints = []
for i in range(count):
    waypoint = sst.getWaypoint(i)
    waypoints.append([waypoint.t, waypoint.x, waypoint.y, waypoint.z, waypoint.vx, waypoint.vy, waypoint.vz, waypoint.ax, waypoint.ay, waypoint.az])

waypoints = np.array(waypoints)
print(waypoints)

ax = visual.visualizeWorld(show=False)
ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'b')
ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'ro')
ax.set_xlabel("-x")
ax.set_ylabel("y")
ax.set_zlabel("-z")
plt.show()
