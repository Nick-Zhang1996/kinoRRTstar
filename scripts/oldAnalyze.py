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

# kinoRRT stuff
world = World(-0.1,3, -3,3, -3,0)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)

# create a small window
obs1 = Box(1,1.3, -3,-2, -3,0)
obs2 = Box(1,1.3, -2,1, -0.7,0)
#obs3 = Box(1,1.3, -2,1, -3,-1.5)
obs4 = Box(1,1.3, 1,3, -3,0)
start_node = Node(0,0,-0.5)
goal_node = Node(2.5,1,-0.5)



#obstacles = [obs1, obs2, obs3, obs4]
obstacles = [obs1, obs2, obs4]
for obstacle in obstacles:
    world.addObstacle(obstacle)
    visual.addObstacle(obstacle)


# TODO check start and goal are in world bobundary
print_info("initializing kinoRRT*")
rrt = KinoRrtStar(world, start_node, goal_node, 600, 10)

try:
    print_info("running kinoRRT*")
    t0 = time()
    rrt.run()
    elapsed = time()-t0
except (RuntimeError):
    print_error("rrt error")

print_info("rrt search finished in " + str(elapsed) +"sec")

waypoint_n = rrt.prepareSolution()

waypoints = []
for i in range(waypoint_n):
    p = rrt.getNextWaypoint()
    assert (p.valid)
    waypoints.append( (p.t,p.x,p.y,p.z) )

# list of (t,x,y,z) 
waypoints = np.array(waypoints)

ax = visual.visualizeWorld(show=False)
ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'b')
ax.scatter(-waypoints[:,1], waypoints[:,2], -waypoints[:,3], 'ro')
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
