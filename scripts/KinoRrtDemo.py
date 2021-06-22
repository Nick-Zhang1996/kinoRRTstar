import os
import sys
import numpy as np
from time import time
import matplotlib.pyplot as plt

base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src/')
sys.path.append(base_dir)
from kinoRRT import *

from WorldVisualization import WorldVisualization


world = World(-10,10,-5,5,-10,0)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)

obs1 = Box(6-10,0-5,0-10,6+2-10,5-5,10-10)
obs2 = Box(6-10,5-5,0-10,6+2-10,5+5-5,5-10)
obs3 = Box(12-10,5-5,0-10,12+2-10,5+5-5,10-10)
obs4 = Box(12-10,0-5,5-10,12+2-10,5-5,5+5-10)
'''
world = World(0,20,0,10,0,10)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)

obs1 = Box(6,0,0,6+2,5,10)
obs2 = Box(6,5,0,6+2,5+5,5)
obs3 = Box(12,5,0,12+2,5+5,10)
obs4 = Box(12,0,5,12+2,5,5+5)
'''

obstacles = [obs1, obs2, obs3, obs4]
for obstacle in obstacles:
    world.addObstacle(obstacle)
    visual.addObstacle(obstacle)


start_node = Node(2-10, 2-5, 2-10)
goal_node = Node(18-10, 8-5, 8-10)

# TODO check start and goal are in world bobundary
rrt = KinoRrtStar(world, start_node, goal_node, 600, 10)

try:
    t0 = time()
    rrt.run()
    elapsed = time()-t0
except (RuntimeError):
    print("rrt error")

print("rrt search finished in " + str(elapsed) +"sec")

waypoint_n = rrt.prepareSolution()

waypoints = []
for i in range(waypoint_n):
    p = rrt.getNextWaypoint()
    assert (p.valid)
    waypoints.append( (p.t,p.x,p.y,p.z) )

# list of (t,x,y,z) 
waypoints = np.array(waypoints)
#print("waypoints")
#print(waypoints)
print("start")
print(waypoints[0])
print("goal")
print(waypoints[-1])

# find max speed and scale time respectively
diff = np.diff(waypoints, axis=0)
v = (diff[:,1]**2+diff[:,1]**2+diff[:,1]**2)**0.5 / diff[:,0]
max_v = 0.5
scale = np.max(v) /  max_v
waypoints[:,0] *= scale

diff = np.diff(waypoints, axis=0)
new_v = (diff[:,1]**2+diff[:,1]**2+diff[:,1]**2)**0.5 / diff[:,0]

ax = visual.visualizeWorld(show=False)
ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'b')
ax.scatter(-waypoints[:,1], waypoints[:,2], -waypoints[:,3], 'ro')
ax.set_xlabel("-x")
ax.set_ylabel("y")
ax.set_zlabel("-z")
plt.show()

