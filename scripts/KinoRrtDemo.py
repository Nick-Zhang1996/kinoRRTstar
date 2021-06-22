import os
import sys
import numpy as np
from time import time
import matplotlib.pyplot as plt

base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src/')
sys.path.append(base_dir)
from kinoRRT import *

from WorldVisualization import WorldVisualization


world = World(-0.1,3, -3,3, -3,0)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)

# create a small window
obs1 = Box(1,1.3, -3,-2, -3,0)
obs2 = Box(1,1.3, -2,1, -0.7,0)
#obs3 = Box(1,1.3, -2,1, -3,-1.5)
obs4 = Box(1,1.3, 1,3, -3,0)
start_node = Node(0,0,-0.5)
goal_node = Node(2.5,2,-0.5)

'''
world = World(-10,10,-5,5,-10,0)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)

obs1 = Box(-4,-2, -5,0, -10,0)
obs2 = Box(-4,-2, 0,5, -10,-5)
obs3 = Box(2,4, 0,5, -10,0)
obs4 = Box(2,4, -5,0, -5,0)
start_node = Node(2-10, 2-5, 2-10)
goal_node = Node(18-10, 8-5, 8-10)
'''


#obstacles = [obs1, obs2, obs3, obs4]
obstacles = [obs1, obs2, obs4]
for obstacle in obstacles:
    world.addObstacle(obstacle)
    visual.addObstacle(obstacle)



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
print("waypoints")
print(waypoints)
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

