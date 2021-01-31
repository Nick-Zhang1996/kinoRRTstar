# kinodynamic RRT star
import numpy as np
from math import log

from World import World
from DoubleIntegrator import DoubleIntegrator
from Tree import Tree,n

class kinoRRT:
    #
    # nodes: nodes to sample
    # if set to 0, quit after first solution is found
    def __init__(self, world, start_pos, end_pos, nodes=0):
        self.segment_len = 4.0

        #pos, vel, cost
        start_state = np.hstack([start_pos,[0,0,0],0]).astype(np.float)
        self.end_pos = np.array(end_pos,dtype=np.float)

        # init tree
        if (nodes == 0):
            self.tree = Tree(start_state,max_node=6000)
            self.max_nodes = 6000
        else:
            self.tree = Tree(start_state,max_node=nodes)
            self.max_nodes = nodes

        self.world = world
        self.di = DoubleIntegrator()

        # for determining best neighbour distance
        self.gamma = 40
        # problem dimension, 3D world: 3
        self.dim = 3
        self.min_neighbour_radius = 3.0

        # solutions contains information about nodes that can reach end_pos
        # each solution is an entry in self.solutions in for [node_id, total_cost_to_endnode]
        self.solutions = []
        self.found_solution = False
        return

    # run kinodynamic RRT for specified start and finish position
    # velocity at initial_pos is 0
    # velocity at final_pos is free
    # cost function is encoded in DoubleIntegrator class in the polynomials
    # and cannot be changed
    def run(self):
        for i in range(self.max_nodes):
            # randomly sample until we add a valid node to tree
            while (not self.sampleSpace()):
                pass
            print(self.tree.n_nodes)


            # have we found a solution?
            if (self.found_solution and self.max_nodes==0):
                break

    # from original implementation
    def calcNeighbourDist(self):
        nun = self.tree.n_nodes
        ner = self.gamma * ( log(nun + 1.0) / nun )**(1.0 / self.dim);
        return min(ner,self.min_neighbour_radius)


    # sample space, if result is successful return True
    # if fail return False
    def sampleSpace(self):

        # randomly sample in world
        new_pos = np.random.rand(3) * np.array(self.world.size)
        # find nearest point
        closest_id, neighbour_id_vec = self.tree.findClosest(new_pos, self.calcNeighbourDist())
        # if too far move to that neighbour
        dist = self.tree.dist2(new_pos,closest_id)
        if ( dist > self.segment_len):
            new_pos = self.tree.getNodePos(closest_id) + self.segment_len * self.tree.vecFromNode(new_pos,closest_id)

            neighbour_id_vec = self.tree.getNeighbour(new_pos, self.calcNeighbourDist())
            dist = self.segment_len

        # collision?
        if (not self.world.checkNoCollision(new_pos.reshape(1,3))):
            # better luck next time
            return False

        # now we have a collision free point
        # and we know its nearest neighbour

        # can we go from neighbour to new node?
        init_state = self.tree.getNodeState(closest_id)
        tf = self.di.DI3d_timeFreeVel(*init_state, *new_pos)
        state_vec = self.di.DI3d_stateFreeVel(tf,*init_state, *new_pos)
        # no need to check start and finish
        interior_state_vec = state_vec[1:-1]
        if (not self.world.checkNoCollision(interior_state_vec)):
            return False

        # no collision, add this sampled point to our tree
        # search neighbour for best parent
        best_parent = closest_id
        min_cost = self.di.DI3d_costFreeVel(tf, *init_state, *new_pos) + self.tree.getNodeCost(closest_id)
        # pos: pos
        # state: pos,vel
        # full_state: pos,vel,cost
        new_full_state = np.hstack([state_vec[-1],min_cost])
        for neighbour in neighbour_id_vec:
            init_state = self.tree.getNodeState(closest_id)

            tf = self.di.DI3d_timeFreeVel(*init_state, *new_pos)
            state_vec = self.di.DI3d_stateFreeVel(tf,*init_state, *new_pos)
            # no need to check start and finish
            interior_state_vec = state_vec[1:-1]
            if (not self.world.checkNoCollision(interior_state_vec)):
                continue
            cost = self.di.DI3d_costFreeVel(tf, *init_state, *new_pos) + self.tree.getNodeCost(neighbour)
            if (cost < min_cost):
                min_cost = cost
                best_parent = neighbour
                new_full_state = np.hstack([state_vec[-1],min_cost])

        # now add new_node to tree
        new_node_id = self.tree.addNode(new_full_state, best_parent)
        new_state = new_full_state[:-1]

        # is the new node close to end_state?
        if (self.tree.dist2(end_pos,new_node_id) < self.segment_len):
            # can we reach end_state?
            tf = self.di.DI3d_time(*new_state, *end_pos)
            state_vec = self.di.DI3d_state(tf,*new_state, *end_pos)
            interior_state_vec = state_vec[1:-1]
            if (self.world.checkNoCollision(interior_state_vec)):
                # we've found a path
                self.tree.setEndState(new_node_id,1)
                cost = self.di.DI3d_costFreeVel(tf, *new_state, *self.end_pos) + self.tree.getNodeCost(new_node_id)
                self.solutions.append([new_node_id,cost])
                self.found_solution = True
                # no need to rewire, new_node should remain a leaf node
                return True

        # rewire
        # can new node be a better parent for its neighbours? 
        for neighbour in neighbour_id_vec:
            neighbour_state = self.tree.getNodeState(neighbour)
            tf = self.di.DI3d_time(*new_state, *neighbour_state)
            state_vec = self.di.DI3d_state(tf,*new_state, *neighbour_state)
            interior_state_vec = state_vec[1:-1]
            if (not self.world.checkNoCollision(interior_state_vec)):
                continue
            cost = self.di.DI3d_cost(tf, *new_state, *neighbour_state) + self.tree.getNodeCost(new_node_id)
            if (cost < self.tree.getNodeCost(neighbour)):
                d_cost = cost - self.tree.getNodeCost(neighbour)

                self.tree.transferChild(new_node_id,neighbour)
                self.tree.spreadCostDelta(neighbour,d_cost)
        return True


    # extract solutions from built tree
    # plot
    #   environment with obstacles
    #   sampled states
    #   trajectory through sampled states
    # if multiple solutions are found show each sequentially
    def showResult(self):
        return


if __name__=="__main__":
    world = World(size=(20,10,10))
    obs = [[6,8],[0,5],[0,10]]
    world.addObstacle(obs)
    obs = [[6,8],[5,10],[0,6]]
    world.addObstacle(obs)

    start_pos = (0,0,0)
    end_pos = (10,3,5)
    # keep searching until finding a solution
    rrt = kinoRRT(world,start_pos, end_pos)
    rrt.run()
    print(rrt.solutions)
