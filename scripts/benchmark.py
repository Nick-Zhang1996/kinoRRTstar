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
import pickle

def makeWorld():
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
    return world, visual, start_node, goal_node

def benchmarkKinoRRT(world, start_node, goal_node, duration):
    print_info("initializing kinoRRT*")
    interior_points = 50
    rrt = KinoRrtStar(world, start_node, goal_node, interior_points)

    try:
        duration = 30
        print_info("running kinoRRT* for %.2f seconds"%(duration))
        t0 = time()
        rrt.runWithTimeLimit(duration)
        elapsed = time()-t0
    except (RuntimeError):
        print_error("rrt error")

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
    '''
    print("key_waypoints")
    print(key_waypoints)
    print("start")
    print(key_waypoints[0])
    print("goal")
    print(key_waypoints[-1])
    '''

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
    '''
    print_info("total time : %.1f sec "%(traj_t))
    print_info("max speed : %.1f m/s "%(max_speed))
    print_info("max acc : %.1f m/s "%(max_acc))
    '''

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

    # progress statistics
    nodes = rrt.getNodeCountHistPy()
    cost = rrt.getMinCostHistPy()
    sols = rrt.getSolutionCountHistPy()
    return np.array(nodes), np.array(cost), np.array(sols)

def benchmarkSST(world, start_node, goal_node, duration):
    print("solving with SST, time limit = %.2f sec"%(duration))
    sst = OmplBenchmark(world, start_node, goal_node)
    sst.solveIncrementally(duration, 1.0)
    count = sst.getWaypointCount()
    waypoints = []
    for i in range(count):
        waypoint = sst.getWaypoint(i)
        waypoints.append([waypoint.t, waypoint.x, waypoint.y, waypoint.z, waypoint.vx, waypoint.vy, waypoint.vz, waypoint.ax, waypoint.ay, waypoint.az])

    waypoints = np.array(waypoints)
    #print(waypoints)

    ax = visual.visualizeWorld(show=False)
    ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'b')
    ax.plot(-waypoints[:,1], waypoints[:,2], -waypoints[:,3],'ro')
    ax.set_xlabel("-x")
    ax.set_ylabel("y")
    ax.set_zlabel("-z")
    plt.show()
    # progress statistics
    nodes = sst.getNodeCountHistPy()
    cost = sst.getMinCostHistPy()
    sols = sst.getSolutionCountHistPy()
    return np.array(nodes), np.array(cost), np.array(sols)

if __name__ == "__main__":
    duration = 60.0
    world, visual, start_node, goal_node = makeWorld()
    rrt_nodes, rrt_cost, rrt_sols = benchmarkKinoRRT(world, start_node, goal_node, duration)
    sst_nodes, sst_cost, sst_sols = benchmarkSST(world, start_node, goal_node, duration)
    pass
    with open("log.txt", 'wb') as f:
        rrt = [rrt_nodes, rrt_cost, rrt_sols]
        sst = [sst_nodes, sst_cost, sst_sols]
        data = [rrt, sst]
        pickle.dump(data,f)

