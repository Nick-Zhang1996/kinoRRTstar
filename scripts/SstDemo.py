
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

np.set_printoptions(precision=4)

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
duration = 10.0
print("solving with SST, time limit = %.2f sec"%(duration))
sst = OmplBenchmark(world, start_node, goal_node)
sst.solveIncrementally(duration, 10.0)

# progress statistics
'''
print(sst.getNodeCountHistPy())
print(sst.getMinCostHistPy())
print(sst.getSolutionCountHistPy())
'''


count = sst.getWaypointCount()
waypoints = []
for i in range(count):
    waypoint = sst.getWaypoint(i)
    #waypoints.append([waypoint.t, waypoint.x, waypoint.y, waypoint.z, waypoint.vx, waypoint.vy, waypoint.vz, waypoint.ax, waypoint.ay, waypoint.az])
    waypoints.append([waypoint.t, waypoint.x, waypoint.y, waypoint.z, waypoint.vx, waypoint.vy, waypoint.vz])

waypoints = np.array(waypoints)
#print(waypoints)

# verify cost
interior_points = 50
rrt = KinoRrtStar(world, start_node, goal_node, interior_points)
'''
cost_oc = 0.0
cost_direct = 0.0
for i in range(waypoints.shape[0]-1):
    p0 = waypoints[i]
    p1 = waypoints[i+1]
    t_s = (p1[1] - p0[1]) / (p1[4] + p0[4]) * 2;
    #print("time %.3f"%(t_s))
    cost_oc += rrt.cost(t_s, p0[1], p0[2], p0[3],p0[4], p0[5], p0[6],p1[1], p1[2], p1[3],p1[4], p1[5], p1[6])
    acc = waypoints[i+1] - waypoints[i]
    acc = acc[4:] / t_s
    print(waypoints[i])
    print(acc)
    print(t_s)
    segment_cost = np.linalg.norm(acc)**2* t_s + t_s
    print_ok(segment_cost)
    cost_direct += segment_cost
print("[python] total cost_oc = %.3f"%(cost_oc))
print("[python] total cost_direct = %.3f"%(cost_direct))
'''


ax = visual.visualizeWorld(show=False)
ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'b')
ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'ro')
ax.set_xlabel("-x")
ax.set_ylabel("y")
ax.set_zlabel("-z")
plt.show()
