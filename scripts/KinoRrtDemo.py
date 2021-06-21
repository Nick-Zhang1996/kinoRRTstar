import os
import sys
import numpy as np
from time import time
import matplotlib.pyplot as plt

base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src/')
sys.path.append(base_dir)
from kinoRRT import *

from WorldVisualization import WorldVisualization


world = World(20,10,10)
visual = WorldVisualization(size=(world.x,world.y,world.z))

obs1 = Box(6,0,0,6+2,5,10)
obs2 = Box(6,5,0,6+2,5+5,5)
obs3 = Box(12,5,0,12+2,5+5,10)
obs4 = Box(12,0,5,12+2,5,5+5)

obstacles = [obs1, obs2, obs3, obs4]
for obstacle in obstacles:
    world.addObstacle(obstacle)
    visual.addObstacle(obstacle)


start_node = Node(2, 2, 2)
goal_node = Node(18, 8, 8)

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

ax = visual.visualizeWorld(show=False)
ax.plot(waypoints[:,1], waypoints[:,2], waypoints[:,2],'b')
ax.scatter(waypoints[:,1], waypoints[:,2], waypoints[:,2], 'ro')
plt.show()

