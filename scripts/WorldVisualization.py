# this class handles 
# creation of a world with boundaries and obstacles
# collision checking
# visualization
import numpy as np
import matplotlib.pyplot as plt
from itertools import product, combinations
from mpl_toolkits.mplot3d import Axes3D
from kinoRRT import World

class WorldVisualization:
    # size: ((x_l,x_h),(y_l,y_h),(z_l,z_h)) dimension of world
    def __init__(self,dim):
        self.dim = dim = np.array(dim,dtype=np.float)
        size = np.array(dim[:,1]-dim[:,0]).flatten()
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
    def addObstacle(self,obs):
        obstacle_size = ((obs.x_l, obs.x_h),(obs.y_l, obs.y_h),(obs.z_l, obs.z_h))
        obstacle = np.array(obstacle_size,dtype=np.float)
        assert obstacle.shape == (3,2)
        assert np.all((obstacle[:,1]-obstacle[:,0])>0.0)
        self.obstacles.append(obstacle)
        return


    def visualizeWorld(self,show=True):
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
                #ax.plot_surface(X, Y, Z, rstride=1, cstride=1, **kwargs)
                ax.plot_wireframe(-X, Y, -Z, rstride=1, cstride=1,color='tan')

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        for obstacle in self.obstacles:
            pos = (obstacle[0,0], obstacle[1,0], obstacle[2,0])
            size = (obstacle[0,1]-obstacle[0,0], obstacle[1,1]-obstacle[1,0], obstacle[2,1]-obstacle[2,0])
            plotCubeAt(pos,size,ax,color='tan')


        # force plt to plot entire world
        pos = self.dim[:,0]
        size = self.size
        X, Y, Z = cuboid_data( pos, size )
        #ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1,color='w')
        ax.scatter(-X, Y, -Z,'k')

        if show:
            plt.show()
        else:
            return ax

