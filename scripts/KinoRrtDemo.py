import os
import sys
import numpy as np
from time import time
import matplotlib.pyplot as plt

base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src/')
sys.path.append(base_dir)
from kinoRRT import *

from WorldVisualization import WorldVisualization
from common import *
from window import createWindow


# Original big test world
world = World(-10,10,-5,5,-10,0)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)

obs1 = Box(-4,-2, -5,0, -10,0)
obs2 = Box(-4,-2, 0,5, -10,-5)
obs3 = Box(2,4, 0,5, -10,0)
obs4 = Box(2,4, -5,0, -5,0)
start_node = Node(2-10, 2-5, 2-10)
goal_node = Node(18-10, 8-5, 8-10)
obstacles = [obs1, obs2, obs3, obs4]
'''

# successful experiment 1, short trajectory
world = World(-0.1,3, -3,3, -3,0)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)

# create a small window
obs1 = Box(1,1.3, -3,-2, -3,0)
obs2 = Box(1,1.3, -2,1, -0.7,0)
obs3 = Box(1,1.3, -2,1, -3,-1.5)
obs4 = Box(1,1.3, 1,3, -3,0)
start_node = Node(0,0,-0.5)
goal_node = Node(2.5,1,-0.5)
obstacles = [obs1,obs2,obs3,obs4]
'''

'''
# newly created window obstacle
start_node = Node(-1.8, 0.6, -0.6)
goal_node = Node(1.7,0,-2)

#world = World(-6,6,-2.5,2.5,-3,0)
world = World(-2,2,-0.7,1,-2.5,0)
dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
visual = WorldVisualization(dim)
obstacles = []
obstacles += createWindow((-1,0,0), 0.5, True, world)
obstacles += createWindow((0.5,-1,0), 0.5, False, world)
'''
for obstacle in obstacles:
    world.addObstacle(obstacle)
    visual.addObstacle(obstacle)

visual.visualizeWorld(show=True)

# TODO check start and goal are in world bobundary
print_info("initializing kinoRRT*")
interior_points = 50
rrt = KinoRrtStar(world, start_node, goal_node, 4000, interior_points)

try:
    print_info("running kinoRRT*")
    t0 = time()
    rrt.run()
    elapsed = time()-t0
except (RuntimeError):
    print_error("rrt error")

print_info("rrt search finished in " + str(elapsed) +"sec")

key_waypoint_n = rrt.prepareSolution()
if (key_waypoint_n == 0):
    print_error("RRT failed to find a solution")
    exit(1)
key_waypoints = []

#for i in range(key_waypoint_n):
p = rrt.getNextWaypoint()
key_waypoints.append( (p.t,p.x,p.y,p.z) )
while (p.valid):
    p = rrt.getNextWaypoint()
    key_waypoints.append( (p.t,p.x,p.y,p.z) )

# list of (t,x,y,z) 
key_waypoints = np.array(key_waypoints)
print("key_waypoints")
print(key_waypoints)
print("start")
print(key_waypoints[0])
print("goal")
print(key_waypoints[-1])

# get continuous waypoint
#print_ok("-------- continuous traj -----")
traj_t = rrt.getTrajectoryTime()
tt = np.linspace(0,traj_t)
waypoints = []
for this_t in tt:
    waypoint = rrt.getTrajectory(this_t)
    if (not waypoint.valid):
        print_error("waypoint error")
        break;
    waypoints.append([waypoint.t, waypoint.x, waypoint.y, waypoint.z, waypoint.vx, waypoint.vy, waypoint.vz, waypoint.ax, waypoint.ay, waypoint.az])
waypoints = np.array(waypoints)
max_speed = (np.max(waypoints[:,4]**2 + waypoints[:,5]**2 + waypoints[:,6]**2))**0.5
max_acc = (np.max(waypoints[:,7]**2 + waypoints[:,8]**2 + waypoints[:,9]**2))**0.5
print_info("total time : %.1f sec "%(traj_t))
print_info("max speed : %.1f m/s "%(max_speed))
print_info("max acc : %.1f m/s "%(max_acc))

# find max speed and scale time respectively
diff = np.diff(waypoints, axis=0)
v = (diff[:,1]**2+diff[:,1]**2+diff[:,1]**2)**0.5 / diff[:,0]
max_v = 0.5
scale = np.max(v) /  max_v
waypoints[:,0] *= scale

diff = np.diff(waypoints, axis=0)
new_v = (diff[:,1]**2+diff[:,1]**2+diff[:,1]**2)**0.5 / diff[:,0]

ax = visual.visualizeWorld(show=False)
#ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'b')
ax.scatter(-key_waypoints[:,1], key_waypoints[:,2], -key_waypoints[:,3], 'ro')
ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'r')
ax.set_xlabel("-x")
ax.set_ylabel("y")
ax.set_zlabel("-z")
plt.show()

