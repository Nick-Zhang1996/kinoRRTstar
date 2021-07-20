#create box obstacles for window obstacles made from insulation
from WorldVisualization import WorldVisualization
from common import *
import matplotlib.pyplot as plt

import os
import sys
base_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../src/')
sys.path.append(base_dir)
from kinoRRT import *

# create a list of box obstacles to make a window ,facing x direction
# the window will also block to the full extent of the world dimension
# leaving only a small opening
# origin: min x,y,-z coordinate of the window
# thickness: in meters
# low: bool, whether opening is low or high
def createWindow(origin, thickness, low, world):
    f2m = 0.3048
    in2m = 0.0254
    # create the actual obstacle 
    if (low):
        # lowest
        obs1 = Box(origin[0], origin[0]+thickness, origin[1], origin[1]+4*f2m, origin[2]-16*in2m, origin[2])
        # window left
        obs2 = Box(origin[0], origin[0]+thickness, origin[1], origin[1]+16*in2m, origin[2]-32*in2m, origin[2]-16*in2m)
        # window right
        obs3 = Box(origin[0], origin[0]+thickness, origin[1]+4*f2m-16*in2m, origin[1]+4*f2m, origin[2]-32*in2m, origin[2]-16*in2m)
        # upper part
        obs4 = Box(origin[0], origin[0]+thickness, origin[1], origin[1]+4*f2m, origin[2]-8*f2m, origin[2]-32*in2m)
    else:
        # lower part
        obs1 = Box(origin[0], origin[0]+thickness, origin[1], origin[1]+4*f2m, origin[2]-8*f2m+32*in2m, origin[2])
        # window left
        obs2 = Box(origin[0], origin[0]+thickness, origin[1], origin[1]+16*in2m, origin[2]-8*f2m+16*in2m, origin[2]-8*f2m+32*in2m)
        # window right
        obs3 = Box(origin[0], origin[0]+thickness, origin[1]+4*f2m-16*in2m, origin[1]+4*f2m,origin[2]-8*f2m+16*in2m, origin[2]-8*f2m+32*in2m)
        # upper
        obs4 = Box(origin[0], origin[0]+thickness, origin[1], origin[1]+4*f2m, origin[2]-8*f2m, origin[2]-8*f2m+16*in2m)
    obstacles = [obs1, obs2, obs3, obs4]

    # create the virtual obstacle to block the remaining open space
    x_l = origin[0]
    x_h = origin[0]+thickness
    y_l = origin[1]
    y_h = origin[1]+4*f2m
    z_l = origin[2]-8*f2m
    z_h = origin[2]

    '''
    if (x_l > world.x_l):
        obs = Box(world.x_l, x_l, world.y_l, world.y_h, world.z_l, world.z_h)
        obstacles.append(obs)
    if (x_h < world.x_h):
        obs = Box(x_h, world.x_h, world.y_l, world.y_h, world.z_l, world.z_h)
        obstacles.append(obs)
    '''
    if (y_l > world.y_l):
        obs = Box(x_l, x_h, world.y_l, y_l, z_l, z_h)
        obstacles.append(obs)
    if (y_h < world.y_h):
        obs = Box(x_l, x_h, y_h, world.y_h, z_l, z_h)
        obstacles.append(obs)
    if (z_l > world.z_l):
        obs = Box(x_l, x_h, world.y_l, world.y_h, world.z_l, z_l)
        obstacles.append(obs)
    if (z_h < world.z_h):
        obs = Box(x_l, x_h, world.y_l, world.y_h, z_h, world.z_h)
        obstacles.append(obs)
    return obstacles


if __name__=="__main__":
    world = World(-6,6,-2.5,2.5,-3,0)
    dim = ((world.x_l, world.x_h), (world.y_l, world.y_h),(world.z_l, world.z_h))
    visual = WorldVisualization(dim)
    obstacles = createWindow((-2,0,0), 1, True, world)
    obstacles += createWindow((2,-1,0), 1, False, world)
    for obs in obstacles:
        visual.addObstacle(obs)

    ax = visual.visualizeWorld(show=False)
    plt.show()
