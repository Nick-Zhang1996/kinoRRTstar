# this class handles 
# creation of a world with boundaries and obstacles
# collision checking
# visualization
import numpy as np
import matplotlib.pyplot as plt
from itertools import product, combinations
from mpl_toolkits.mplot3d import Axes3D

class World:
    # size: (x,y,z) size of world
    def __init__(self,size):
        size = np.array(size,dtype=np.float).flatten()
        assert size.shape[0] == 3
        assert np.all(size>0)
        self.size = size

        # obstacle
        # each obstacle is rectangular box represented by a 3*2 array
        # ((x_low_bound, x_high_bound), (y_l, y_h), (z_l,z_h))
        self.obstacles = []

    # add obstacle to world
    # each obstacle is rectangular box represented by a 3*2 array
    # obstacle_size : ((x_low_bound, x_high_bound), (y_l, y_h), (z_l,z_h))
    def addObstacle(self,obstacle_size):
        obstacle = np.array(obstacle_size,dtype=np.float)
        assert obstacle.shape == (3,2)
        assert np.all((obstacle[:,1]-obstacle[:,0])>0.0)
        self.obstacles.append(obstacle)
        return

    # given a list-like of (x,y,z), check if any of the coordinates is in collision
    # or out side of boundary
    # coord: np.array of size (n,3+), n being number of coordinates in the batch, only first three states are used (x,y,z,...)
    # NOTE does not check data type and format
    # return True if no cllision, false if collide
    def checkNoCollision(self,coord):
        # boundary
        if (np.any(coord<0.0)):
            return False
        if (np.any(self.size[0] < coord[:,0]) or
            np.any(self.size[1] < coord[:,1]) or 
            np.any(self.size[2] < coord[:,2])):
            return False

        # check obstacle
        for obstacle in self.obstacles:
            x_cond = np.logical_and(coord[:,0] > obstacle[0,0], coord[:,0] < obstacle[0,1])
            y_cond = np.logical_and(coord[:,1] > obstacle[1,0], coord[:,1] < obstacle[1,1])
            z_cond = np.logical_and(coord[:,2] > obstacle[2,0], coord[:,2] < obstacle[2,1])
            cond = np.logical_and(x_cond,y_cond)
            cond = np.logical_and(cond,z_cond)
            if np.any(cond):
                return False
                break

        return True

    def visualizeWorld(self):
        def cuboid_data(o, size=(1,1,1)):
            # code taken from
            # https://stackoverflow.com/a/35978146/4124317
            # suppose axis direction: x: to left; y: to inside; z: to upper
            # get the length, width, and height
            l, w, h = size
            x = [[o[0], o[0] + l, o[0] + l, o[0], o[0]],
                 [o[0], o[0] + l, o[0] + l, o[0], o[0]],
                 [o[0], o[0] + l, o[0] + l, o[0], o[0]],
                 [o[0], o[0] + l, o[0] + l, o[0], o[0]]]
            y = [[o[1], o[1], o[1] + w, o[1] + w, o[1]],
                 [o[1], o[1], o[1] + w, o[1] + w, o[1]],
                 [o[1], o[1], o[1], o[1], o[1]],
                 [o[1] + w, o[1] + w, o[1] + w, o[1] + w, o[1] + w]]
            z = [[o[2], o[2], o[2], o[2], o[2]],
                 [o[2] + h, o[2] + h, o[2] + h, o[2] + h, o[2] + h],
                 [o[2], o[2], o[2] + h, o[2] + h, o[2]],
                 [o[2], o[2], o[2] + h, o[2] + h, o[2]]]
            return np.array(x), np.array(y), np.array(z)

        def plotCubeAt(pos=(0,0,0), size=(1,1,1), ax=None,**kwargs):
            # Plotting a cube element at position pos
            if ax !=None:
                X, Y, Z = cuboid_data( pos, size )
                ax.plot_surface(X, Y, Z, rstride=1, cstride=1, **kwargs)

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        for obstacle in self.obstacles:
            pos = (obstacle[0,0], obstacle[1,0], obstacle[2,0])
            size = (obstacle[0,1]-obstacle[0,0], obstacle[1,1]-obstacle[1,0], obstacle[2,1]-obstacle[2,0])
            plotCubeAt(pos,size,ax,color='b')


        # force plt to plot entire world
        pos = (0,0,0)
        size = self.size
        X, Y, Z = cuboid_data( pos, size )
        #ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1,color='w')
        ax.scatter(X, Y, Z,'k')

        plt.show()



if __name__=="__main__":
    test_world = World(size=(20,10,10))
    obs = [[6,8],[0,5],[0,10]]
    test_world.addObstacle(obs)
    obs = [[6,8],[5,10],[0,6]]
    test_world.addObstacle(obs)

    p1 = np.array([[7,4,6]]) # collide
    p2 = np.array([[5,4,6]])
    p3 = np.array([[5,8,6]])
    pnts = np.vstack([p1,p2,p3])
    # expect False
    print(test_world.checkNoCollision(pnts))
    pnts = np.vstack([p1,p2])
    # expect False
    print(test_world.checkNoCollision(pnts))
    pnts = np.vstack([p3,p2])
    # expect True
    print(test_world.checkNoCollision(pnts))

    # test plotting
    test_world.visualizeWorld()






        

